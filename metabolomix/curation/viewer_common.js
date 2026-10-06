const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const el = (tag, attrs = {}, text) => { const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text; return e; };

function preparerFor(c) {
  const ratio = window.devicePixelRatio || 1;
  const prepare = () => { const w = c.clientWidth, h = c.clientHeight; c.width = w * ratio; c.height = h * ratio;
    const ctx = c.getContext("2d"); ctx.scale(ratio, ratio); return {ctx, w, h}; };
  prepare.canvas = c;
  return prepare;
}

function canvasFor(card) {
  const c = el("canvas"); card.appendChild(c);
  return preparerFor(c);
}

const fixed = (v, digits) => typeof v === "number" && isFinite(v) ? v.toFixed(digits) : "–";

function drawEic(prepare, spot) {
  const {ctx, w, h} = prepare(); const samples = spot.eic.samples; if (!samples.length) return;
  const xs = samples.flatMap(s => s.points.map(p => p[0])); const ys = samples.flatMap(s => s.points.map(p => p[1]));
  if (!xs.length) return;
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y1 = Math.max(...ys) || 1;
  const px = x => 4 + (x - x0) / ((x1 - x0) || 1) * (w - 8), py = y => h - 12 - y / y1 * (h - 18);
  const rep = samples.find(s => s.representative) || samples[0];
  ctx.fillStyle = css("--line"); ctx.globalAlpha = 0.5;
  ctx.fillRect(px(rep.left), 4, px(rep.right) - px(rep.left), h - 16); ctx.globalAlpha = 1;
  for (const s of samples) {
    ctx.beginPath();
    ctx.strokeStyle = s.partner ? css("--ref") : s.representative ? css("--meas") : css("--muted");
    ctx.lineWidth = s.representative || s.partner ? 2 : 1;
    ctx.setLineDash(s.partner ? [2, 2] : s.detected ? [] : [4, 3]);
    s.points.forEach((p, i) => i ? ctx.lineTo(px(p[0]), py(p[1])) : ctx.moveTo(px(p[0]), py(p[1])));
    ctx.stroke();
  }
  ctx.setLineDash([]); ctx.fillStyle = css("--muted"); ctx.font = "10px system-ui";
  ctx.fillText(x0.toFixed(2) + "–" + x1.toFixed(2) + " min", 4, h - 2);
}

// --- mirror labels (pure) ---
// library_plot_mirror の `_draw_labels`（auto 方式）と同じ貪欲法: 高さの降順に走査し、既に置いた
// ラベルと横・縦の両方が近いものを飛ばす。文字列は m/z（小数 4 桁）、片側 25 本が上限。
// 側（実測・参照）ごとに別々に呼ぶ——側をまたいだ衝突は見ない（上流も別 Annotator）。
const MIRROR_MAX_LABELS_PER_SIDE = 25;

// `bounds`（{xmin, xmax}、任意）を渡すとラベルの中心をその内側に寄せる（端のピークのラベルが
// 縦軸の目盛りやキャンバスの外にはみ出さないように）。重なりは寄せた後の位置で判定する。
function pickMirrorLabels(points, toPixel, measure, bounds) {
  const placed = [];
  for (const [mz, height] of [...points].sort((a, b) => b[1] - a[1])) {
    if (placed.length >= MIRROR_MAX_LABELS_PER_SIDE) break;
    const text = mz.toFixed(4), box = measure(text), p = {...toPixel(mz, height)};
    if (bounds) p.x = Math.min(Math.max(p.x, bounds.xmin + box.w / 2), bounds.xmax - box.w / 2);
    if (placed.some(o => (box.w + o.box.w) / 2 > Math.abs(p.x - o.x)
                      && (box.h + o.box.h) / 2 > Math.abs(p.y - o.y))) continue;
    placed.push({text, mz, x: p.x, y: p.y, box});
  }
  return placed;
}
// --- end mirror labels ---

function drawMirror(prepare, spot) {
  const {ctx, w, h} = prepare(); const m = spot.mirror;
  if (!m) { ctx.fillStyle = css("--muted"); ctx.font = "12px system-ui";
    ctx.fillText(spot.mirror_empty_text || (spot.match && !spot.match.has_msms ? "MS/MS なし" : "参照スペクトルなし"), 8, h / 2); return; }
  const all = m.measured.concat(m.reference).map(p => p[0]);
  const left = 26, headroom = 12;       // 縦軸目盛りの幅、ラベルの逃げ
  const x0 = Math.min(...all) - 5, x1 = Math.max(...all) + 5, mid = h / 2, amp = mid - headroom;
  const px = x => left + (x - x0) / ((x1 - x0) || 1) * (w - left - 4);
  const norm = pts => { const top = Math.max(...pts.map(p => p[1])) || 1; return pts.map(p => [p[0], p[1] / top]); };
  const measured = norm(m.measured), reference = norm(m.reference);
  const matched = new Set(m.matched_mz.map(v => v.toFixed(2)));
  for (const [x, y] of measured) { ctx.strokeStyle = css("--meas"); ctx.beginPath();
    ctx.moveTo(px(x), mid); ctx.lineTo(px(x), mid - y * amp); ctx.stroke(); }
  for (const [x, y] of reference) { ctx.strokeStyle = matched.has(x.toFixed(2)) ? css("--match") : css("--ref");
    ctx.beginPath(); ctx.moveTo(px(x), mid); ctx.lineTo(px(x), mid + y * amp); ctx.stroke(); }
  ctx.strokeStyle = css("--line"); ctx.beginPath(); ctx.moveTo(left, mid); ctx.lineTo(w, mid);
  ctx.moveTo(left, mid - amp); ctx.lineTo(left, mid + amp); ctx.stroke();
  // 縦軸: 相対強度 %（実測は上、参照は下。それぞれ自分の最大値を 100 とする）
  ctx.fillStyle = css("--muted"); ctx.font = "9px system-ui"; ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const pct of [0, 50, 100]) for (const sign of pct ? [1, -1] : [1]) {
    const y = mid - sign * pct / 100 * amp;
    ctx.fillText(String(pct), left - 4, y); ctx.beginPath(); ctx.moveTo(left - 2, y); ctx.lineTo(left, y); ctx.stroke();
  }
  ctx.textAlign = "center";
  const measure = text => ({w: ctx.measureText(text).width, h: 9});
  for (const [pts, sign] of [[measured, 1], [reference, -1]]) {
    ctx.textBaseline = sign > 0 ? "bottom" : "top";
    const toPixel = (mz, y) => ({x: px(mz), y: mid - sign * y * amp});
    for (const label of pickMirrorLabels(pts, toPixel, measure, {xmin: left, xmax: w}))
      ctx.fillText(label.text, label.x, label.y - sign * 2);
  }
}
