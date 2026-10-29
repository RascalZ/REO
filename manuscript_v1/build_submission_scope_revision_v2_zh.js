// Local direct translation renderer. Reads only the translated manuscript and original figure PNGs.
const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, ImageRun,
  WidthType, BorderStyle, AlignmentType, HeadingLevel, Footer, Header, PageNumber,
} = require('docx');
const root = path.resolve(__dirname, '..');
const out = path.join(__dirname, 'submission_final_v2_zh');
const stem = 'Manuscript_BMC_Genomics_submission_v2_zh';
const lines = fs.readFileSync(path.join(out, stem + '.md'), 'utf8').split(/\r?\n/);
const body = [];
const font = { ascii: 'Times New Roman', hAnsi: 'Times New Roman', eastAsia: '宋体', cs: 'Times New Roman' };
const headFont = { ...font, eastAsia: '黑体' };
const WIDTH = 11906, HEIGHT = 16838, MARGIN = 1247, BODY_WIDTH = WIDTH - 2 * MARGIN;
const figureFiles = [null, 'Figure1_study_design.png', 'Figure2_development_primary.png',
  'Figure3_representation_and_zero_solutions.png', 'Figure4_loto_and_external.png'];
const captions = new Map(lines.filter(l => /^图[1-4]\. /.test(l)).map(l => [Number(l[1]), l]));
const anchors = [null, '图1概述了队列、嵌套分区、模型臂及评价边界。',
  '图2展示配对的各药物及等权药物估计。', '图3B展示全开发集拟合的分子系数状态。',
  '图4B展示外部各药物和宏差值。'];
const emitted = new Set();
let section = '', tableCount = 0;

function runs(text, extra = {}) {
  const result = [], re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g;
  let last = 0;
  for (const m of text.matchAll(re)) {
    if (m.index > last) result.push(new TextRun({ text: text.slice(last, m.index), font, size: 22, ...extra }));
    const bold = m[0].startsWith('**');
    result.push(new TextRun({ text: m[0].slice(bold ? 2 : 1, bold ? -2 : -1), font, size: 22,
      ...(bold ? { bold: true } : { italics: true }), ...extra }));
    last = m.index + m[0].length;
  }
  if (last < text.length) result.push(new TextRun({ text: text.slice(last), font, size: 22, ...extra }));
  return result.length ? result : [new TextRun({ text, font, size: 22, ...extra })];
}
function para(text, options = {}) {
  return new Paragraph({ children: runs(text), spacing: { line: 360, after: 100 },
    widowControl: true, ...options });
}
function figure(n) {
  const bytes = fs.readFileSync(path.join(root, 'manuscript_v1/submission_final_v2/Figures', figureFiles[n]));
  const nw = bytes.readUInt32BE(16), nh = bytes.readUInt32BE(20), w = 620;
  body.push(new Paragraph({
    children: [new ImageRun({ type: 'png', data: bytes, transformation: { width: w, height: Math.round(w * nh / nw) },
      altText: { title: `图${n}`, description: captions.get(n), name: figureFiles[n] } })],
    alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60, line: 240 }, keepNext: true,
  }));
  const caption = captions.get(n);
  const stop = caption.indexOf('。');
  const formatted = stop >= 0 ? `**${caption.slice(0, stop + 1)}**${caption.slice(stop + 1)}` : caption;
  body.push(new Paragraph({ children: runs(formatted, { size: 19 }),
    spacing: { line: 270, after: 160 }, keepLines: true }));
  emitted.add(n);
}
function table(rows) {
  tableCount++;
  const widths = tableCount === 1 ? [1630, 850, 1850, 850, BODY_WIDTH - 5180]
    : [1650, 820, 2700, 850, BODY_WIDTH - 6020];
  const border = { style: BorderStyle.SINGLE, size: 2, color: '808080' };
  return new Table({ width: { size: BODY_WIDTH, type: WidthType.DXA }, columnWidths: widths,
    borders: { top: border, bottom: border, left: border, right: border, insideHorizontal: border, insideVertical: border },
    rows: rows.map((row, ri) => new TableRow({ tableHeader: ri === 0, cantSplit: true,
      children: row.map((value, ci) => new TableCell({ width: { size: widths[ci], type: WidthType.DXA },
        margins: { top: 60, bottom: 60, left: 65, right: 65 },
        children: [new Paragraph({ children: runs(value.replace(/ to /g, ' 至 '), { size: 17, bold: ri === 0 }),
          spacing: { line: 240, before: 10, after: 10 }, keepNext: true })] })) })) });
}

body.push(new Paragraph({ children: [new TextRun({ text: '正式英文修订v2的中文直接翻译', font: headFont, size: 30, bold: true })],
  alignment: AlignmentType.CENTER, spacing: { after: 120 }, keepNext: true }));
body.push(new Paragraph({ children: [new TextRun({ text: '译稿日期：2026年10月29日。以下按原稿逐段翻译，不改写或增删研究内容。中文版独立编页码，反馈请注明“中文PDF第X页”。药物名、作者姓名、参考文献题录及图内标签保留原文，图注和表头译为中文。', font, size: 18, color: '666666' })],
  spacing: { line: 270, after: 220 }, keepLines: true }));
for (let i = 0; i < lines.length;) {
  const line = lines[i].trimEnd();
  if (!line.trim()) { i++; continue; }
  if (line.startsWith('# ')) {
    body.push(new Paragraph({ children: runs(line.slice(2), { size: 22, bold: true }),
      alignment: AlignmentType.CENTER, spacing: { after: 150 }, keepNext: true }));
    i++; continue;
  }
  if (line.startsWith('## ')) {
    section = line.slice(3);
    if (section !== '标题页' && section !== '图注') body.push(new Paragraph({
      children: runs(section, { size: 26, bold: true, font: headFont }), heading: HeadingLevel.HEADING_1,
      spacing: { before: 240, after: 120, line: 360 }, keepNext: true }));
    i++; continue;
  }
  if (line.startsWith('### ')) {
    body.push(new Paragraph({ children: runs(line.slice(4), { size: 23, bold: true, font: headFont }),
      heading: HeadingLevel.HEADING_2, spacing: { before: 160, after: 80, line: 340 }, keepNext: true }));
    i++; continue;
  }
  if (section === '图注') { i++; continue; }
  if (line.startsWith('|')) {
    const rows = [];
    while (i < lines.length && lines[i].startsWith('|')) {
      const row = lines[i].trim().replace(/^\|/, '').replace(/\|$/, '').split('|').map(x => x.trim());
      if (!row.every(c => /^:?-{2,}:?$/.test(c))) rows.push(row);
      i++;
    }
    body.push(table(rows)); continue;
  }
  const text = line.replace(/`/g, '');
  let p = para(text);
  if (text.startsWith('**作者：**')) p = new Paragraph({ children: [
    new TextRun({ text: '作者：', font, size: 22, bold: true }),
    new TextRun({ text: text.slice('**作者：**'.length), font, size: 22 }),
  ], spacing: { line: 360, after: 100 }, keepNext: true });
  if (text.startsWith('**题目：**')) p = new Paragraph({ children: runs(text, { size: 27, bold: true }),
    spacing: { line: 400, after: 160 }, keepNext: true });
  if (/^表[12]\. /.test(text)) {
    const stop = text.indexOf('。');
    body.push(new Paragraph({ children: runs(`**${text.slice(0, stop + 1)}**`, { size: 20 }),
      spacing: { line: 280, before: 120, after: 80 }, keepNext: true, keepLines: true }));
    // Preserve the full original legend following the table, while keeping the title before it.
    const legend = text.slice(stop + 1);
    let j = i + 1;
    while (j < lines.length && !lines[j].trim()) j++;
    const rows = [];
    while (j < lines.length && lines[j].startsWith('|')) {
      const row = lines[j].trim().replace(/^\|/, '').replace(/\|$/, '').split('|').map(x => x.trim());
      if (!row.every(c => /^:?-{2,}:?$/.test(c))) rows.push(row);
      j++;
    }
    if (!rows.length) throw new Error('Missing table rows');
    body.push(table(rows));
    body.push(new Paragraph({ children: runs(legend, { size: 19 }), spacing: { line: 270, before: 80, after: 100 }, keepLines: true }));
    i = j; continue;
  }
  if (text.startsWith('*注：*')) p = new Paragraph({ children: runs(text, { size: 18 }),
    spacing: { line: 270, after: 160 }, keepLines: true });
  if (section === '参考文献') p = new Paragraph({ children: runs(text, { size: 18 }),
    spacing: { line: 260, after: 80 }, indent: { left: 280, hanging: 280 }, keepLines: true });
  if (section === '补充文件') p = new Paragraph({ children: runs(text, { size: 20 }),
    spacing: { line: 300, after: 100 }, keepLines: true });
  body.push(p);
  for (let n = 1; n <= 4; n++) if (!emitted.has(n) && text.includes(anchors[n])) figure(n);
  i++;
}
if (tableCount !== 2 || emitted.size !== 4) throw new Error(`Expected2tables/4figures;got${tableCount}/${emitted.size}`);
const doc = new Document({ creator: 'REO local manuscript translation', title: '正式英文修订v2的中文直译稿',
  description: 'Direct reading translation only; English manuscript remains authoritative and unchanged.',
  styles: { default: { document: { run: { font, size: 22 }, paragraph: { spacing: { line: 360, after: 100 } } } } },
  sections: [{ properties: { page: { size: { width: WIDTH, height: HEIGHT },
    margin: { top: 1150, right: MARGIN, bottom: 1150, left: MARGIN, header: 550, footer: 550 } } },
    headers: { default: new Header({ children: [new Paragraph({ alignment: AlignmentType.RIGHT,
      children: [new TextRun({ text: '中文通读 · 源：正式英文修订v2 · 2026-10-29', font, size: 16, color: '777777' })] })] }) },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ children: ['第 ', PageNumber.CURRENT, ' 页 / 共 ', PageNumber.TOTAL_PAGES, ' 页'], font, size: 18 })] })] }) },
    children: body }] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync(path.join(out, stem + '.docx'), b); console.log('Built Chinese reading Word with2tables/4originalPNGs.'); })
  .catch(e => { console.error(e); process.exit(1); });
