"""Paso A5: versión 3D del conducto mamario (cilindro, no círculo).

Por qué importa el 3D y no es un capricho estético: en 2D toda célula del tumor está
cerca del borde, así que cualquier terapia la alcanza. En 3D una masa tumoral es una
ESFERA, y su interior queda lejos de todo. La relación superficie/volumen cae con el
radio (∝ 1/r), de modo que cuanto más crece, menor fracción de sus células es accesible.
Esa es la razón principal por la que terapias que funcionan en placa fracasan en un tumor.

Reutiliza la misma biología y las mismas políticas que la versión 2D: aquí solo cambia
la geometría, la difusión y el vecindario (26 vecinos en vez de 8).

    python tejido3d.py --dias 400 --terapia guardiana+memoria --dia-terapia 200
"""
from __future__ import annotations

import argparse
import random
import time
from dataclasses import replace

import numpy as np

from agentes import Agente, Analizar, Atacar, Mover, Observacion, Percepcion, Reparar
from biologia import ESCENARIOS, GENES, PASOS_POR_DIA, Parametros
from politicas import TERAPIAS
from tejido import Celula, Simulacion, LUMINAL, MIOEPITELIAL, FIBROBLASTO

LUZ, EPITELIO, MIO, ESTROMA = range(4)
NOMBRE_ZONA = ("luz", "epitelio", "mioepitelio", "estroma")
VECINOS = tuple((dz, dy, dx)
                for dz in (-1, 0, 1) for dy in (-1, 0, 1) for dx in (-1, 0, 1)
                if dz or dy or dx)
MAX_AGENTES = 1200
VIDA_NECROSIS = 40


def _signo(v: float) -> int:
    return (v > 0) - (v < 0)


class Simulacion3D:
    """Conducto mamario como cilindro hueco a lo largo del eje Z."""

    def __init__(self, escenario: str = "luminal", semilla: int | None = None,
                 params: Parametros | None = None, lado: int = 54, largo: int = 40,
                 radio: float = 13.0):
        esc = ESCENARIOS[escenario]
        self.escenario = escenario
        self.p = replace(params or Parametros(), estrogeno=esc.estrogeno)
        self.rng = random.Random(semilla)
        self.t = 0
        self._siguiente_id = 0
        self.forma = (largo, lado, lado)
        self.celulas: dict[tuple, Celula] = {}
        self.agentes: list[Agente] = []
        self.stats = dict(apoptosis_p53=0, anoikis=0, necrosis=0, muertes_inmunes=0,
                          muertes_terapia=0, colateral=0, reparaciones=0,
                          reparaciones_inutiles=0, roturas_membrana=0, vasos_nuevos=0,
                          desdiferenciaciones=0, guardianas=0, fagocitadas=0,
                          fagocitadas_sanas=0, dosis=0, madres_sanas_tocadas=0, propagaciones=0,
                          pasos_bloqueados=0)
        self._construir(esc, radio)
        for _ in range(40):
            self._difundir()

    # ------------------------------------------------------------------ montaje
    def _nuevo_id(self) -> int:
        self._siguiente_id += 1
        return self._siguiente_id

    def _construir(self, esc, radio: float):
        largo, lado, _ = self.forma
        self.centro = (lado - 1) / 2
        zz, yy, xx = np.mgrid[0:largo, 0:lado, 0:lado]
        self.dist = np.hypot(yy - self.centro, xx - self.centro)   # distancia al eje del cilindro
        self.zona = np.full(self.forma, ESTROMA, np.uint8)
        self.zona[self.dist < radio] = MIO
        self.zona[self.dist < radio - 1.5] = EPITELIO
        self.zona[self.dist < radio - 3.2] = LUZ

        azar = self.rng
        self.vasos = np.zeros(self.forma, bool)
        periductal = (self.dist >= radio + 2) & (self.dist < radio + 3.5)
        lejos = self.dist >= radio + 3.5
        ruido = np.array(azar.random()) * 0 + np.random.default_rng(
            azar.randrange(1 << 30)).random(self.forma)
        self.vasos[(periductal & (ruido < 0.05)) | (lejos & (ruido < 0.012))] = True

        for pos in zip(*np.nonzero(self.zona == EPITELIO)):
            self.celulas[tuple(map(int, pos))] = Celula(
                self._nuevo_id(), LUMINAL, dict(esc.germinal), azar.random() < esc.fraccion_er,
                madre=azar.random() < self.p.fraccion_madre_normal)
        for pos in zip(*np.nonzero(self.zona == MIO)):
            self.celulas[tuple(map(int, pos))] = Celula(self._nuevo_id(), MIOEPITELIAL)
        for pos in zip(*np.nonzero((self.zona == ESTROMA) & ~self.vasos & (ruido < 0.04))):
            self.celulas[tuple(map(int, pos))] = Celula(self._nuevo_id(), FIBROBLASTO)

        if esc.semilla:
            madres = [q for q, c in self.celulas.items() if c.tipo == LUMINAL and c.madre]
            c = self.celulas[azar.choice(madres)]
            c.genoma.update(esc.semilla)
            c.neo = sum(azar.random() < self.p.prob_neoantigeno for _ in esc.semilla)
            c.actualizar_fenotipo()

        self.O2 = np.full(self.forma, 0.8)
        self.Q = np.zeros(self.forma)
        self.P = np.zeros(self.forma)   # presión intersticial: empuja hacia fuera de la masa
        self.factor_presion = 1.0
        self.cd47_hasta = -1
        self.vasos_tumorales = np.zeros(self.forma, bool)
        self.modo_vasos = "normal"

    # ------------------------------------------------------------------ bucle
    @property
    def dia(self) -> float:
        return self.t / PASOS_POR_DIA

    def _vecinos(self, pos):
        largo, lado, _ = self.forma
        z, y, x = pos
        for dz, dy, dx in VECINOS:
            nz, ny, nx = z + dz, y + dy, x + dx
            if 0 <= nz < largo and 0 <= ny < lado and 0 <= nx < lado:
                yield nz, ny, nx

    def paso(self):
        self.t += 1
        orden = list(self.celulas.items())
        self.rng.shuffle(orden)
        for pos, c in orden:
            if self.celulas.get(pos) is c:
                self._vivir(pos, c)
        self._difundir()
        self._mover_agentes()

    def _vivir(self, pos, c: Celula):
        rng, p = self.rng, self.p
        if c.tipo in (FIBROBLASTO, MIOEPITELIAL):
            return
        if c.necrosis:
            c.necrosis += 1
            if c.necrosis > VIDA_NECROSIS:
                del self.celulas[pos]
            return
        o2 = self.O2[pos]
        # MUERTE POR O2 GRADUAL (variante B). Misma regla que tejido.py (ver allí la
        # justificación completa: Durand y Sham IJROBP 1998;42:711, y Papandreou et al.
        # Cancer Res 2005;65:3171). Hasta el 2026-10-29 este módulo seguía con la regla
        # vieja —sólo morían las de O2 < o2_necrosis, con p = 0,08— que hacía a las
        # hipóxicas casi inmortales. Las dos mallas tienen que envejecer igual o ninguna
        # medida 3D es comparable con la 2D.
        if o2 < p.o2_necrosis:
            p_necrosis = 0.08
        elif o2 < p.o2_hipoxia:
            p_necrosis = 0.08 * (p.o2_hipoxia - o2) / (p.o2_hipoxia - p.o2_necrosis)
        else:
            p_necrosis = 0.0
        if p_necrosis and rng.random() < p_necrosis:
            c.necrosis = 1
            self.stats["necrosis"] += 1
            return
        if c.n_drivers == 0 and rng.random() < p.prob_muerte_normal * (0.1 if c.madre else 1):
            del self.celulas[pos]
            return
        if c.n_drivers and not c.madre:
            if c.divisiones == 0 and rng.random() < p.prob_muerte_diferenciada:
                del self.celulas[pos]
                return
            sasp = any((v := self.celulas.get(q)) is not None and v.senescente
                       for q in self._vecinos(pos))
            if rng.random() < (p.plasticidad * (3 if o2 < p.o2_hipoxia else 1)
                               + (p.prob_sasp_plasticidad if sasp else 0)):
                c.madre, c.divisiones = True, None
                self.stats["desdiferenciaciones"] += 1

        z = self.zona[pos]
        fuera = (z == LUZ and not c.anoikis_res) or (z in (MIO, ESTROMA) and not c.invasiva)
        if fuera and rng.random() < p.prob_anoikis:
            del self.celulas[pos]
            self.stats["anoikis"] += 1
            return
        if c.invasiva and z != ESTROMA and rng.random() < p.prob_invasion * 8:
            q = next(iter(self._vecinos(pos)))
            v = self.celulas.get(q)
            if v is not None and v.tipo == MIOEPITELIAL:
                del self.celulas[q]
                self.stats["roturas_membrana"] += 1
        if c.guardiana and rng.random() < p.prob_fagocitosis * (
                p.coste_lectura_interna if c.guardiana in Simulacion.PUERTAS_INTERNAS else 1.0):
            self._fagocitar(pos, c)
            if self.celulas.get(pos) is not c:
                return
        if c.senescente:
            return
        pdiv = p.prob_division * c.factor * ((1 + p.estrogeno) if c.er else 1)
        if o2 < p.o2_hipoxia:
            pdiv *= 0.15
        if c.madre and c.n_drivers:
            pdiv *= p.factor_quiescencia_madre
        if rng.random() < pdiv:
            self._dividir(pos, c)

    def _permitido(self, c: Celula, origen, destino) -> bool:
        zd = self.zona[destino]
        if self.vasos[destino]:
            return False
        if zd == EPITELIO:
            return True
        if zd == LUZ:
            return c.anoikis_res
        if zd == MIO:
            return c.invasiva
        return c.invasiva and self.zona[origen] in (MIO, ESTROMA)

    def _dividir(self, pos, c: Celula):
        if c.divisiones == 0:
            return
        huecos = [q for q in self._vecinos(pos)
                  if q not in self.celulas and self._permitido(c, pos, q)]
        if not huecos and c.senal_prolif:
            rivales = [q for q in self._vecinos(pos)
                       if (v := self.celulas.get(q)) is not None and v.tipo == LUMINAL
                       and v.n_drivers == 0 and not v.necrosis and self._permitido(c, pos, q)]
            if rivales and self.rng.random() < 0.3:
                del self.celulas[rivales[0]]
                huecos = rivales[:1]
        if not huecos:
            return
        q = self.rng.choice(huecos)
        hija = Celula(self._nuevo_id(), LUMINAL, dict(c.genoma), c.er, c.mutaciones)
        hija.reparada, hija.neo = c.reparada, c.neo
        hija.cd47 = c.cd47 if self.t < self.cd47_hasta else True
        pa = self.p.prob_aneuploidia if not c.p53_ok else self.p.prob_aneuploidia_p53ok
        hija.aneuploide = c.aneuploide or self.rng.random() < pa
        hija.telomerasa_ok = c.telomerasa_ok
        inmortal = (c.inmortal or c.madre) and c.telomerasa_ok
        if c.madre and c.n_drivers:
            hija.madre = self.rng.random() < self.p.prob_division_simetrica
            hija.divisiones = None if (hija.madre and inmortal) else self.p.divisiones_obrera
            hija.guardiana = c.guardiana if hija.madre else None
        elif c.madre:
            hija.madre = self.rng.random() < self.p.fraccion_madre_normal
            hija.guardiana = c.guardiana if hija.madre else None
        elif inmortal:
            hija.divisiones = c.divisiones = None
        else:
            if c.divisiones is None:
                c.divisiones = self.p.divisiones_obrera
            c.divisiones -= 1
            hija.divisiones = c.divisiones
        for cel, sitio in ((c, pos), (hija, q)):
            golpes = self._mutar(cel)
            dano = golpes + (1 if cel.hrd and self.rng.random() < self.p.prob_dano_hrd else 0)
            if dano and cel.p53_ok and \
                    self.rng.random() < 1 - (1 - self.p.prob_apoptosis_p53) ** dano:
                self.stats["apoptosis_p53"] += 1
                if cel is c:
                    del self.celulas[pos]
                continue
            if cel is hija:
                self.celulas[q] = hija

    def _mutar(self, c: Celula) -> int:
        mult = self.p.tasa_mutacion * (self.p.factor_hrd if c.hrd else 1) * (1 if c.p53_ok else 3)
        golpes = 0
        for nombre, gen in GENES.items():
            if self.rng.random() < mult * gen.peso:
                h = c.genoma.get(nombre, 0)
                if h < (1 if gen.tipo == "oncogen" else 2):
                    c.genoma[nombre] = h + 1
                    golpes += 1
                    c.neo += self.rng.random() < self.p.prob_neoantigeno
        lam = 0.3 * (10 if c.hrd else 1)
        nuevas = int(lam) + (self.rng.random() < lam - int(lam))
        c.mutaciones += nuevas
        c.neo += sum(self.rng.random() < self.p.prob_neoantigeno / 15 for _ in range(nuevas))
        if golpes:
            c.actualizar_fenotipo()
            if c.n_drivers and not c.madre and c.divisiones is None and not c.inmortal:
                c.divisiones = self.p.divisiones_obrera
        return golpes

    # ------------------------------------------------------------------ campos
    def _difundir(self, iteraciones: int = 4):
        consumo = np.zeros(self.forma)
        fuente = np.zeros(self.forma)
        for pos, c in self.celulas.items():
            if c.necrosis:
                continue
            consumo[pos] = 0.016 if c.n_drivers else 0.004
            if c.n_drivers:
                fuente[pos] = 1.0
        # Presión intersticial (ver tejido.py): la generan el tumor y los fibroblastos,
        # los linfáticos sanos la drenan fuera. En 3D es donde de verdad importa: la
        # meseta de presión cubre todo el interior de la esfera.
        if self.p.presion_intersticial:
            pres = fuente.copy()
            for pos, c in self.celulas.items():
                if c.tipo == FIBROBLASTO:
                    pres[pos] = 0.4
            pres *= self.p.presion_intersticial * self.factor_presion
            drenaje = np.where(pres > 0, 0.02, 0.25)
            for _ in range(iteraciones):
                lap = -6 * self.P
                for eje in (0, 1, 2):
                    lap = lap + np.roll(self.P, 1, eje) + np.roll(self.P, -1, eje)
                self.P += 0.14 * lap + 0.08 * pres - drenaje * self.P
                np.clip(self.P, 0.0, 1.0, out=self.P)
        for campo, sumidero, aporte in ((self.O2, consumo, None), (self.Q, 0.08, fuente)):
            for _ in range(iteraciones):
                lap = -6 * campo
                for eje in (0, 1, 2):
                    lap = lap + np.roll(campo, 1, eje) + np.roll(campo, -1, eje)
                campo += 0.14 * lap - sumidero * campo
                if aporte is None:
                    campo[self.vasos] = 1.0
                else:
                    campo += 0.05 * aporte
                np.clip(campo, 0.0, 1.0, out=campo)

    # ------------------------------------------------------------------ agentes
    def _observar(self, c: Celula, pos, dz, dy, dx) -> Observacion:
        hipoxica = self.O2[pos] < self.p.o2_hipoxia
        neo = c.neo if (c.mhc and c.n_drivers > 0) else 0
        tejido = ("luminal", "mioepitelial", "fibroblasto")[c.tipo]
        return Observacion(dy, dx, c.id, tejido, c.mhc, neo, c.pdl1, c.her2,
                           c.senescente or c.n_drivers >= 2 or (hipoxica and c.n_drivers > 0),
                           self.t - c.ultima_division < 8, hipoxica, c.madre, c.senescente,
                           c.cd47, Simulacion._expone_comeme(c, hipoxica),
                           c.madre and not self._en_nicho(pos))

    def _percibir(self, ag: Agente) -> Percepcion:
        z, y, x = ag.pos
        self._dz = {}
        vecinas = []
        for dz in range(-1, 2):
            for dy in range(-2, 3):
                for dx in range(-2, 3):
                    c = self.celulas.get((z + dz, y + dy, x + dx))
                    if c is not None and not c.necrosis:
                        o = self._observar(c, (z + dz, y + dy, x + dx), dz, dy, dx)
                        self._dz[(o.dy, o.dx)] = dz
                        vecinas.append(o)
        mejor, quimio = self.Q[z, y, x] + 1e-4, None
        for q in self._vecinos(ag.pos):
            if self.Q[q] > mejor:
                mejor, quimio = self.Q[q], (q[1] - y, q[2] - x)
        return Percepcion(vecinas, quimio, (_signo(self.centro - y), _signo(self.centro - x)),
                          NOMBRE_ZONA[self.zona[z, y, x]], float(self.O2[z, y, x]), self.dia,
                          self.rng)

    def _mover_agentes(self):
        largo, lado, _ = self.forma
        vivos = []
        for ag in self.agentes:
            ag.vida -= 1
            if ag.vida <= 0 or ag.agotamiento >= 1:
                continue
            vivos.append(ag)
            p = self._percibir(ag)
            accion = ag.politica(ag, p)
            if isinstance(accion, Mover):
                dz = self.rng.randint(-1, 1)
                destino = (min(largo - 1, max(0, ag.pos[0] + dz)),
                           min(lado - 1, max(0, ag.pos[1] + _signo(accion.dy))),
                           min(lado - 1, max(0, ag.pos[2] + _signo(accion.dx))))
                if Simulacion._frena_la_presion(self, ag.pos, destino):
                    self.stats["pasos_bloqueados"] += 1
                else:
                    ag.pos = destino
                continue
            if not isinstance(accion, (Atacar, Analizar, Reparar)):
                continue
            dz = self._dz.get((accion.dy, accion.dx), 0)
            q = (ag.pos[0] + dz, ag.pos[1] + accion.dy, ag.pos[2] + accion.dx)
            c = self.celulas.get(q)
            if c is None or c.necrosis:
                continue
            if isinstance(accion, Atacar):
                self._atacar(ag, c, q, vivos)
            elif isinstance(accion, Analizar):
                ag.memoria.setdefault("genomas", {})[c.id] = {g: h for g, h in c.genoma.items() if h}
            else:
                self._reparar(ag, c, q, accion.gen)
        self.agentes = vivos

    def _reconoce(self, ag: Agente, c: Celula, q) -> str | None:
        r = ag.receptores
        if "INESPECIFICO" in r:
            return "inespecifico"
        if c.tipo != LUMINAL:
            return None
        o = self._observar(c, q, 0, 0, 0)
        if "HER2" in r and c.her2:
            return "car"
        if "TCR" in r and o.neoantigenos > 0:
            return "tcr"
        if "NK" in r and (not c.mhc or o.ligandos_estres):
            return "nk"
        return None

    def _atacar(self, ag: Agente, c: Celula, q, vivos: list):
        via = self._reconoce(ag, c, q)
        if via is None:
            return
        if via in ("tcr", "car") and c.pdl1 and not ag.resistente_pdl1:
            ag.agotamiento += 0.2
            self.stats["frenos_pdl1"] = self.stats.get("frenos_pdl1", 0) + 1
            # DESENGANCHE: el mismo artefacto que se corrigió en tejido.py el
            # 2026-10-28 y que aquí seguía sin corregir hasta el 2026-10-29.
            # agentes.politica_linfocito_t YA consultaba ag.memoria["frenado_por"],
            # pero este módulo no lo escribía nunca: sin memoria, cazar() vuelve en el
            # paso siguiente a la MISMA célula PD-L1+ y el T se queda anclado a ella
            # hasta agotarse. La señal PD-1 inhibe al CTL pero no lo ancla — al
            # contrario, bloquea la señal de parada del TCR y el T sigue moviéndose
            # (Fife, Nat Immunol 2009;10:1185). El freno sigue costando agotamiento y
            # el T puede volver a esa célula pasado el desenganche si sigue ahí.
            ag.memoria.setdefault("frenado_por", {})[c.id] = \
                self.dia + self.p.dias_desenganche_pd1
            return
        if self.rng.random() >= self.p.prob_muerte_ataque:
            return
        del self.celulas[q]
        ag.muertes += 1
        self.stats["muertes_inmunes" if ag.tipo in ("T", "NK") else "muertes_terapia"] += 1
        if c.tipo != LUMINAL or c.n_drivers == 0:
            self.stats["colateral"] += 1
        if ag.propaga:
            vecinas = [w for w in self._vecinos(q)
                       if (u := self.celulas.get(w)) is not None and u.tipo == LUMINAL
                       and u.n_drivers > 0 and not u.necrosis]
            self.rng.shuffle(vecinas)
            for w in vecinas[:int(ag.propaga)]:
                if len(vivos) + 1 >= MAX_AGENTES:
                    break
                vivos.append(replace(ag, pos=w, memoria={}, agotamiento=0.0, muertes=0,
                                     reparaciones=0))
                self.stats["propagaciones"] += 1
        if ag.expansion and len(vivos) + 1 < MAX_AGENTES and self.rng.random() < ag.expansion:
            vivos.append(replace(ag, memoria={}, agotamiento=0.0, muertes=0, reparaciones=0))

    def _reparar(self, ag: Agente, c: Celula, q, gen: str):
        if gen not in ag.vectores or ag.carga <= 0 or c.tipo != LUMINAL:
            return
        ag.carga -= 1
        if not c.genoma.get(gen):
            self.stats["reparaciones_inutiles"] += 1
            return
        if self.rng.random() >= self.p.eficacia_reparacion:
            return
        tenia_p53_roto = not c.p53_ok
        c.genoma[gen] = 0
        c.actualizar_fenotipo()
        ag.reparaciones += 1
        self.stats["reparaciones"] += 1
        if c.n_drivers == 0:
            c.reparada = True
            if c.madre and ag.reprograma:
                c.guardiana = ag.reprograma
                self.stats["guardianas"] += 1
        elif tenia_p53_roto and c.p53_ok and self.rng.random() < 0.5:
            del self.celulas[q]
            self.stats["apoptosis_p53"] += 1

    def _en_nicho(self, pos) -> bool:
        """Capa basal: epitelio en contacto con el mioepitelio (ver tejido.py)."""
        if self.zona[pos] != EPITELIO:
            return False
        return any(self.zona[q] == MIO for q in self._vecinos(pos))

    def _fagocitar(self, pos, c: Celula):
        for q in self._vecinos(pos):
            v = self.celulas.get(q)
            if v is None or v.tipo != LUMINAL or v.necrosis or v.guardiana:
                continue
            o = self._observar(v, q, 0, 0, 0)
            diana = Simulacion._puerta(self, c.guardiana, v, o)
            if not diana:
                continue
            del self.celulas[q]
            self.stats["fagocitadas"] += 1
            if v.n_drivers == 0:
                self.stats["fagocitadas_sanas"] += 1
                self.stats["colateral"] += 1
            if self.rng.random() < 0.5:
                self._dividir(pos, c)
            return

    # ------------------------------------------------------------------ terapias
    def _vaso_aleatorio(self):
        zs, ys, xs = np.nonzero(self.vasos)
        i = self.rng.randrange(len(zs))
        return int(zs[i]), int(ys[i]), int(xs[i])

    def _posiciones_entrega(self, n: int, via: str):
        """De dónde salen las células terapéuticas.

        'vena'          por la sangre: aparecen en un vaso y tienen que abrirse paso (lo de siempre)
        'intratumoral'  aguja guiada por ecografía: aparecen DENTRO de la masa, repartidas
        'aguja'         una sola punción: todas en un punto del tumor (lo realista con una aguja)
        """
        tumor = [q for q, c in self.celulas.items()
                 if c.tipo == LUMINAL and c.n_drivers and not c.necrosis]
        if via == "vena" or not tumor:
            return [self._vaso_aleatorio() for _ in range(n)]
        if via == "aguja":
            centro = self.rng.choice(tumor)
            cerca = sorted(tumor, key=lambda q: sum((a - b) ** 2 for a, b in zip(q, centro)))[:n]
            return [self.rng.choice(cerca) for _ in range(n)]
        return [self.rng.choice(tumor) for _ in range(n)]

    def liberar_terapia(self, nombre: str, cantidad: int | None = None, via: str = "vena"):
        if "+" in nombre:
            for parte in nombre.split("+"):
                self.liberar_terapia(parte, cantidad, via)
            return
        cfg = TERAPIAS[nombre]
        n = cantidad or cfg["n"]
        self.stats["dosis"] += n
        for pos in self._posiciones_entrega(n, via):
            ag = cfg["crear"]()
            ag.pos = pos
            self.agentes.append(ag)

    def tratar_vasos(self, modo: str):
        """Mismos modos que en 2D. En 3D solo cambia el efecto sobre el riego; la caída
        de presión al normalizar es la misma (cerrar las fugas baja la PIF)."""
        self.modo_vasos = modo
        if modo == "normalizar":
            self.factor_presion *= 1 - self.p.caida_presion_normalizar
        # OJO: el 3D todavía NO modela la malla (CAFs, colágeno, LOX) ni, por tanto,
        # `tratar_matriz`. Solo tiene presión y geometría. degradar_matriz aquí se limita
        # a retirar fibroblastos.

    # ------------------------------------------------------------------ sin bisturí
    # Mismas reglas que tejido.py (ablacionar, ctdna, focos_carcinoma), en 3D.
    # Estudio: docs/avances/2026-11-05-tratar-sin-operar.md
    def es_carcinoma(self, pos, c: Celula) -> bool:
        """El mismo criterio que 'peligrosas' en resumen(): in situ, invasivo o reina."""
        return (c.tipo == LUMINAL and not c.necrosis and c.n_drivers > 0
                and (c.madre or c.anoikis_res or self.zona[pos] in (MIO, ESTROMA)))

    def focos_carcinoma(self) -> list[list]:
        """Componentes conexas (26 vecinas) de células de carcinoma."""
        celdas = {q for q, c in self.celulas.items() if self.es_carcinoma(q, c)}
        focos, vistos = [], set()
        for q in celdas:
            if q in vistos:
                continue
            foco, pila = [], [q]
            vistos.add(q)
            while pila:
                z, y, x = pila.pop()
                foco.append((z, y, x))
                for dz, dy, dx in VECINOS:
                    v = (z + dz, y + dy, x + dx)
                    if v in celdas and v not in vistos:
                        vistos.add(v)
                        pila.append(v)
            focos.append(foco)
        return focos

    def ablacionar(self, cobertura: float = 0.95, margen: int = 2, borde: int = 1,
                   supervivencia_borde: float = 0.3, guiada: bool = True,
                   umbral_visible: int = 1) -> dict:
        """Ablación en 3D. guiada=True: bola de radio `margen` alrededor de cada célula de
        carcinoma de los focos visibles (≥ umbral_visible células). guiada=False: una
        esfera alrededor del centroide de todas las células con driver (el «disco»).
        Las luminales destruidas quedan como necrosis; las mioepiteliales y los
        fibroblastos se retiran (en 3D _vivir no limpia su necrosis)."""
        if guiada:
            focos = self.focos_carcinoma()
            dianas = [q for f in focos if len(f) >= umbral_visible for q in f]
            vistos = sum(len(f) >= umbral_visible for f in focos)
            radio_bola = margen
        else:
            dianas = [q for q, c in self.celulas.items()
                      if c.tipo == LUMINAL and c.n_drivers and not c.necrosis]
            vistos = 1 if dianas else 0
        if not dianas:
            return dict(destruidas=0, sanas=0, residuo=0, focos_vistos=0)
        dmin: dict = {}
        if guiada:
            ofs = [(dz, dy, dx) for dz in range(-margen, margen + 1)
                   for dy in range(-margen, margen + 1) for dx in range(-margen, margen + 1)
                   if (dz * dz + dy * dy + dx * dx) ** 0.5 <= radio_bola]
            for (z, y, x) in dianas:
                for dz, dy, dx in ofs:
                    q = (z + dz, y + dy, x + dx)
                    d = (dz * dz + dy * dy + dx * dx) ** 0.5
                    if d < dmin.get(q, 1e9):
                        dmin[q] = d
            limite = margen
        else:
            cz, cy, cx = (sum(q[i] for q in dianas) / len(dianas) for i in range(3))
            dist = sorted(((q[0] - cz) ** 2 + (q[1] - cy) ** 2 + (q[2] - cx) ** 2) ** 0.5
                          for q in dianas)
            limite = dist[min(len(dist) - 1, int(len(dist) * cobertura))] + margen
            for q in self.celulas:
                d = ((q[0] - cz) ** 2 + (q[1] - cy) ** 2 + (q[2] - cx) ** 2) ** 0.5
                if d <= limite:
                    dmin[q] = d
        destruidas = sanas = 0
        for q, d in dmin.items():
            c = self.celulas.get(q)
            if c is None or c.necrosis:
                continue
            if d > limite - borde and self.rng.random() < supervivencia_borde:
                continue
            if c.tipo == LUMINAL and c.n_drivers:
                destruidas += 1
            else:
                sanas += 1
            if c.tipo == LUMINAL:
                c.necrosis = 1
            else:
                del self.celulas[q]
        residuo = sum(1 for q, c in self.celulas.items() if self.es_carcinoma(q, c))
        return dict(destruidas=destruidas, sanas=sanas, residuo=residuo, focos_vistos=vistos)

    def ctdna(self, fragmentos_por_celula: float, peso_necrosis: float = 3.0) -> dict:
        vivas = sum(1 for q, c in self.celulas.items() if self.es_carcinoma(q, c))
        muertas = sum(1 for c in self.celulas.values()
                      if c.tipo == LUMINAL and c.n_drivers and c.necrosis)
        lam = fragmentos_por_celula * (vivas + peso_necrosis * muertas)
        n = int(np.random.default_rng(self.rng.getrandbits(32)).poisson(lam))
        return dict(positiva=n > 0, fragmentos=n, esperado=lam, vivas=vivas)

    def degradar_matriz(self, fraccion: float = 0.25):
        quitados = 0
        for pos, c in list(self.celulas.items()):
            if c.tipo == FIBROBLASTO and self.rng.random() < fraccion:
                del self.celulas[pos]
                quitados += 1
        self.factor_presion *= 1 - self.p.caida_presion_matriz * fraccion
        return quitados

    def danar_telomeros(self):
        muertas = jubiladas = sanas = 0
        for pos, c in list(self.celulas.items()):
            if c.tipo != LUMINAL or not (c.inmortal or c.madre) or not c.telomerasa_ok:
                continue
            if self.rng.random() >= self.p.eficacia_6tiodg:
                continue
            if c.n_drivers == 0 and c.madre:
                sanas += 1
            if c.p53_ok:
                del self.celulas[pos]
                self.stats["apoptosis_p53"] += 1
                muertas += 1
            else:
                c.senescente, c.divisiones = True, 0
                jubiladas += 1
        self.stats["madres_sanas_tocadas"] += sanas
        return muertas, jubiladas, sanas

    # ------------------------------------------------------------------ informe
    def resumen(self) -> dict:
        tumor = peligrosas = reinas = senescentes = normales = 0
        for pos, c in self.celulas.items():
            if c.tipo != LUMINAL or c.necrosis:
                continue
            if c.n_drivers == 0:
                normales += 1
                continue
            tumor += 1
            if c.madre:
                reinas += 1
            if c.senescente or c.divisiones == 0:
                senescentes += 1
            if c.madre or c.anoikis_res or self.zona[pos] in (MIO, ESTROMA):
                peligrosas += 1
        return dict(dia=round(self.dia, 1), tumorales=tumor, peligrosas=peligrosas,
                    reinas=reinas, senescentes=senescentes, normales=normales,
                    agentes=len(self.agentes), **self.stats)

    def fraccion_accesible(self) -> float:
        """Fracción de células tumorales con al menos un hueco adyacente: por ahí, y solo
        por ahí, puede llegar una célula terapéutica. Es la medida clave del 3D."""
        tumor = [q for q, c in self.celulas.items()
                 if c.tipo == LUMINAL and c.n_drivers and not c.necrosis]
        if not tumor:
            return 1.0
        libres = sum(any(w not in self.celulas for w in self._vecinos(q)) for q in tumor)
        return libres / len(tumor)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--dias", type=float, default=400)
    ap.add_argument("--dia-terapia", type=float, default=200)
    ap.add_argument("--terapia", default=None)
    ap.add_argument("--tiodg", action="store_true")
    ap.add_argument("--semillas", type=int, default=1)
    ap.add_argument("--via", default="vena", choices=("vena", "intratumoral", "aguja"))
    ap.add_argument("--dosis", type=int, default=None, help="células por terapia (por defecto, la suya)")
    ap.add_argument("--matriz", type=float, default=0.0, help="fracción de matriz degradada (0-1)")
    ap.add_argument("--neoantigenos", type=float, default=0.1)
    a = ap.parse_args()

    print(f"3D · {a.escenario} · {a.dias:g} días · terapia {a.terapia or 'ninguna'} "
          f"el día {a.dia_terapia:g} · vía {a.via}"
          f"{f' · dosis {a.dosis}' if a.dosis else ''}"
          f"{f' · matriz -{a.matriz:.0%}' if a.matriz else ''}"
          f"{' + 6-tio-dG' if a.tiodg else ''}")
    for semilla in range(a.semillas):
        sim = Simulacion3D(a.escenario, semilla, Parametros(prob_neoantigeno=a.neoantigenos))
        t0 = time.time()
        paso_terapia = int(a.dia_terapia * PASOS_POR_DIA)
        for paso in range(int(a.dias * PASOS_POR_DIA)):
            if paso == paso_terapia:
                if a.matriz:
                    sim.degradar_matriz(a.matriz)
                if a.terapia:
                    sim.liberar_terapia(a.terapia, a.dosis, a.via)
                if a.tiodg:
                    sim.danar_telomeros()
            sim.paso()
        r = sim.resumen()
        print(f"  semilla {semilla}: tumorales {r['tumorales']:6} · peligrosas {r['peligrosas']:6} "
              f"· reinas {r['reinas']:5} · senescentes {r['senescentes']:5} "
              f"· colateral {r['colateral']:4} · accesible {sim.fraccion_accesible():.0%} "
              f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
