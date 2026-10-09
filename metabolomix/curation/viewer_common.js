// ビューアの受け口（curation_review / curation_suggest が立てる 127.0.0.1 の HTTP）。MCP Apps と
// 受け口を立てられなかったときは null（Copy / Send だけ）。spec 2026-10-09。
const SUBMIT_ENDPOINT = /*__SUBMIT_ENDPOINT__*/null;
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
    ctx.fillText(spot.mirror_empty_text || (spot.match && !spot.match.has_msms ? "No MS/MS" : "No reference spectrum"), 8, h / 2); return; }
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

// --- submit client (pure) ---
// ビューアの Submit / Finish（spec 2026-10-09）。受け口の状態: unavailable（受け口なし・未接続）/
// live（ping が通った）/ lost（live の後で 401 か届かない＝時間切れ・サーバ再起動）/ finished。
function countFlags(flags) { const n = {}; for (const f of flags) n[f.flag] = (n[f.flag] || 0) + 1; return n; }
function reviewConfirmText(flags) {
  const n = countFlags(flags);
  return `Record Wrong ${n.wrong || 0} / Suspect ${n.suspect || 0} / Confirmed ${n.confirmed || 0} / clear ${n.clear || 0} and update _tags.xml.\n` +
    "If this project is open in MS-DIAL, close it first.";
}
// assign / redundant は _tags.xml を変えない。clear（元の注釈に戻す）だけが Misannotation と Confirmed を外す。
function suggestConfirmText(flags) {
  const n = countFlags(flags);
  const head = `assign ${n.assign || 0} / redundant ${n.redundant || 0} / clear ${n.clear || 0} を記録します。`;
  return n.clear ? head + "\nclear は _tags.xml の Misannotation と Confirmed を外します。" +
    "MS-DIAL でこのプロジェクトを開いているなら先に閉じてください。" : head;
}
function stateAfter(current, httpStatus) {
  if (current === "finished" || current === "lost") return current;
  if (httpStatus === 401 || httpStatus === 0) return current === "live" ? "lost" : "unavailable";
  if (httpStatus === 200) return "live";
  return current;             // 400 / 409 などは接続の状態を変えない
}
function submitControls(state) {
  const live = state === "live";
  return {submit: live, finish: live, copyPrimary: !live};
}
function tagsResultText(body, lang) {
  body = body || {};
  const t = body.tags_xml || {}, c = t.confirmed || {}, n = body.recorded, len = a => (a || []).length;
  if (t.error) return lang === "ja" ? `${n} 件を記録しました。_tags.xml への反映に失敗: ${t.error}`
                                    : `Recorded ${n}. Updating _tags.xml failed: ${t.error}`;
  const tags = `Misannotation +${len(t.added)} / -${len(t.removed)}, Confirmed +${len(c.added)} / -${len(c.removed)}`;
  return (lang === "ja" ? `${n} 件を記録しました。${tags}。` : `Recorded ${n}. ${tags}. `) + (t.note || "");
}
function lostText(lang) {
  return lang === "ja" ? "送信の受け口がありません（時間切れかサーバの再起動）。下の欄をコピーしてチャットに貼ってください。未送信の選択は残っています。"
    : "The submit endpoint is gone (timed out or the server restarted). Copy the text below and paste it into the chat; your unsent changes are kept.";
}
function finishConfirmText(n, lang) {
  return lang === "ja" ? `${n} 件の未送信の変更があります。終了してよいですか（後から「送信用テキストをコピー」で送れます）。`
    : `${n} unsent change(s). Finish anyway? (You can still send them later with Copy submission text.)`;
}
function finishedText(lang) {
  return lang === "ja" ? "終了しました。以降は「送信用テキストをコピー」で送れます。"
    : "Finished. Use Copy submission text to send anything else.";
}
// --- end submit client ---

const HEARTBEAT_MS = 5 * 60 * 1000;

// 本文は text/plain で送る（事前確認の要らない単純要求。中身は JSON）。届かなければ status 0。
async function postEndpoint(path, extra) {
  if (!SUBMIT_ENDPOINT) return {status: 0, body: null};
  try {
    const res = await fetch(`http://127.0.0.1:${SUBMIT_ENDPOINT.port}${path}`, {method: "POST",
      headers: {"Content-Type": "text/plain;charset=UTF-8"},
      body: JSON.stringify({token: SUBMIT_ENDPOINT.token, ...extra})});
    let body = null; try { body = await res.json(); } catch { /* 本文なし */ }
    return {status: res.status, body};
  } catch { return {status: 0, body: null}; }
}

function startSubmitClient(hooks) {
  let state = "unavailable";
  const update = status => { const next = stateAfter(state, status);
    if (next !== state) { state = next; hooks.onState(state); } };
  const ping = async () => { if (state === "lost" || state === "finished") return;
    update((await postEndpoint("/v1/ping", {})).status); };
  if (SUBMIT_ENDPOINT) {
    ping(); setInterval(ping, HEARTBEAT_MS);
    document.addEventListener("visibilitychange", () => { if (document.visibilityState === "visible") ping(); });
  }
  return {
    state: () => state,
    async submit(payload) { const r = await postEndpoint("/v1/submit", payload); update(r.status); return r; },
    async finish() { const r = await postEndpoint("/v1/finish", {}); state = "finished"; hooks.onState(state); return r; },
  };
}

function showSubmitControls(state) {
  const c = submitControls(state);
  document.getElementById("submit").hidden = !c.submit;
  document.getElementById("finish").hidden = !c.finish;
  const copy = document.getElementById("copy");
  copy.classList.toggle("primary", c.copyPrimary); copy.classList.toggle("quiet", !c.copyPrimary);
}

function openFallback(text, message) {
  const area = document.getElementById("fallback");
  area.hidden = false; area.value = text; area.select();
  document.getElementById("result").textContent = message;   // #status は render() が書き直すので #result に出す
}
