// Lienzo del tejido 2D: pinta la rejilla de estados, el oxígeno, los vasos y los agentes.
// Recibe fotogramas compactos (base64) desde Python con draw(frame).

function decode(b64) {
  const bin = atob(b64);
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

function hexToRgb(hex) {
  if (!hex) return null;
  const v = parseInt(hex.slice(1), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

// mapa de calor del oxígeno: azul (hipoxia) -> transparente (normal)
function o2Color(v) {
  const t = v / 255;
  return [40 + 60 * t, 70 + 140 * t, 220, Math.round(170 * (1 - t))];
}

export default {
  template: `
    <div class="tissue-canvas-wrap">
      <canvas ref="canvas" @click="onClick" @mousemove="onMove" @mouseleave="hover = null"></canvas>
      <div v-if="hover" class="tissue-canvas-hover">{{ hover }}</div>
    </div>`,
  props: {
    palette: Array,
    background: Object,
    vesselColor: String,
    agentStyle: Object,
  },
  data() {
    return { hover: null };
  },
  mounted() {
    this.frame = null;
    this.zone = null;
    this.n = 0;
    this.off = document.createElement("canvas");
    this.rgbPalette = (this.palette || []).map(hexToRgb);
    this.observer = new ResizeObserver(() => this.render());
    this.observer.observe(this.$el);
  },
  unmounted() {
    this.observer.disconnect();
  },
  methods: {
    setZone(n, zoneB64) {
      this.n = n;
      this.zone = decode(zoneB64);
      this.off.width = n;
      this.off.height = n;
    },
    draw(frame) {
      this.frame = {
        states: decode(frame.states),
        vessels: decode(frame.vessels),
        o2: frame.o2 ? decode(frame.o2) : null,
        agents: frame.agents || [],
        selection: frame.selection,
      };
      this.render();
    },
    render() {
      const canvas = this.$refs.canvas;
      if (!canvas || !this.frame || !this.n) return;
      const size = this.$el.clientWidth;
      const dpr = window.devicePixelRatio || 1;
      canvas.width = Math.round(size * dpr);
      canvas.height = Math.round(size * dpr);
      canvas.style.width = size + "px";
      canvas.style.height = size + "px";

      const n = this.n;
      const ctxOff = this.off.getContext("2d");
      const img = ctxOff.createImageData(n, n);
      const lumen = hexToRgb(this.background.lumen);
      const stroma = hexToRgb(this.background.stroma);
      const vessel = hexToRgb(this.vesselColor);
      const { states, vessels, o2 } = this.frame;
      for (let i = 0; i < n * n; i++) {
        let rgb = this.zone[i] <= 1 ? lumen : stroma;
        const s = states[i];
        if (s > 0 && this.rgbPalette[s]) rgb = this.rgbPalette[s];
        if (vessels[i]) rgb = vessel;
        let [r, g, b] = rgb;
        if (o2) {
          const [or, og, ob, oa] = o2Color(o2[i]);
          const a = oa / 255;
          r = r * (1 - a) + or * a;
          g = g * (1 - a) + og * a;
          b = b * (1 - a) + ob * a;
        }
        const k = i * 4;
        img.data[k] = r;
        img.data[k + 1] = g;
        img.data[k + 2] = b;
        img.data[k + 3] = 255;
      }
      ctxOff.putImageData(img, 0, 0);

      const ctx = canvas.getContext("2d");
      ctx.imageSmoothingEnabled = false;
      ctx.drawImage(this.off, 0, 0, canvas.width, canvas.height);

      const cell = canvas.width / n;
      for (const [kind, y, x] of this.frame.agents) {
        const style = this.agentStyle[kind];
        if (!style) continue;
        this.glyph(ctx, style.shape, (x + 0.5) * cell, (y + 0.5) * cell, Math.max(cell * 0.95, 4 * dpr), style.color);
      }
      if (this.frame.selection) {
        const [y, x] = this.frame.selection;
        ctx.strokeStyle = "#ffffff";
        ctx.lineWidth = 1.5 * dpr;
        ctx.beginPath();
        ctx.arc((x + 0.5) * cell, (y + 0.5) * cell, cell * 2.2, 0, 2 * Math.PI);
        ctx.stroke();
      }
    },
    glyph(ctx, shape, cx, cy, s, color) {
      ctx.fillStyle = color;
      ctx.strokeStyle = "rgba(0,0,0,0.65)";
      ctx.lineWidth = Math.max(s * 0.12, 1);
      ctx.beginPath();
      if (shape === "triangle") {
        ctx.moveTo(cx, cy - s * 0.6);
        ctx.lineTo(cx + s * 0.55, cy + s * 0.45);
        ctx.lineTo(cx - s * 0.55, cy + s * 0.45);
      } else if (shape === "diamond") {
        ctx.moveTo(cx, cy - s * 0.6);
        ctx.lineTo(cx + s * 0.5, cy);
        ctx.lineTo(cx, cy + s * 0.6);
        ctx.lineTo(cx - s * 0.5, cy);
      } else if (shape === "star") {
        for (let i = 0; i < 10; i++) {
          const r = i % 2 === 0 ? s * 0.7 : s * 0.3;
          const a = (Math.PI / 5) * i - Math.PI / 2;
          ctx.lineTo(cx + r * Math.cos(a), cy + r * Math.sin(a));
        }
      } else {
        ctx.arc(cx, cy, s * 0.4, 0, 2 * Math.PI);
      }
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    },
    cellAt(ev) {
      if (!this.n) return null;
      const rect = this.$refs.canvas.getBoundingClientRect();
      const x = Math.floor(((ev.clientX - rect.left) / rect.width) * this.n);
      const y = Math.floor(((ev.clientY - rect.top) / rect.height) * this.n);
      if (x < 0 || y < 0 || x >= this.n || y >= this.n) return null;
      return [y, x];
    },
    onClick(ev) {
      const p = this.cellAt(ev);
      if (p) this.$emit("cell", { y: p[0], x: p[1] });
    },
    onMove(ev) {
      const p = this.cellAt(ev);
      this.hover = p ? `(${p[0]}, ${p[1]})` : null;
    },
  },
};
