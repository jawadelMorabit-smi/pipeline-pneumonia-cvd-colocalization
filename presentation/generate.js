#!/usr/bin/env node
'use strict';

/*
 * generate.js — convert any .pptx into a self-contained interactive HTML deck.
 *
 *   node generate.js <input.pptx> [output.html] [options]
 *
 * Options:
 *   --fragments      reveal shapes/bullets step-by-step (opt-in; default: whole slide)
 *   --title "..."    override the document title (default: input file name)
 *   --open           open the result in the default browser when done
 *
 * The output is one HTML file with all text, images and styling inlined — no
 * server, no external assets, works offline and is trivially shareable.
 */

const fs = require('fs');
const path = require('path');
const { execFile } = require('child_process');
const { parsePptx } = require('./lib/pptx-parser');
const { buildHtml } = require('./lib/html-builder');

function parseArgs(argv) {
  const args = { _: [], fragments: false, open: false, title: null };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a === '--fragments') args.fragments = true;
    else if (a === '--open') args.open = true;
    else if (a === '--title') args.title = argv[++i];
    else if (a === '-h' || a === '--help') args.help = true;
    else args._.push(a);
  }
  return args;
}

function usage() {
  console.log(`\nInteractive PPTX converter\n\n  node generate.js <input.pptx> [output.html] [--fragments] [--title "..."] [--open]\n`);
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  if (args.help || args._.length === 0) { usage(); process.exit(args.help ? 0 : 1); }

  const input = path.resolve(args._[0]);
  if (!fs.existsSync(input)) { console.error(`✗ Not found: ${input}`); process.exit(1); }

  const output = args._[1]
    ? path.resolve(args._[1])
    : input.replace(/\.pptx$/i, '') + '.interactive.html';
  const title = args.title || path.basename(input).replace(/\.pptx$/i, '').replace(/[_-]+/g, ' ');

  console.log(`→ Reading  ${path.basename(input)}`);
  const buffer = fs.readFileSync(input);
  const model = await parsePptx(buffer);
  console.log(`  Parsed ${model.slides.length} slides · stage ${Math.round(model.size.w)}×${Math.round(model.size.h)}px`);

  const elementCount = model.slides.reduce((n, s) => n + s.elements.length, 0);
  const withNotes = model.slides.filter((s) => s.notes).length;
  console.log(`  ${elementCount} elements · ${withNotes} slides with speaker notes`);

  const html = buildHtml(model, { title, fragments: args.fragments });
  fs.writeFileSync(output, html, 'utf8');
  const kb = (Buffer.byteLength(html) / 1024).toFixed(0);
  console.log(`✓ Wrote    ${path.basename(output)}  (${kb} KB, self-contained)`);
  console.log(`\n  Open it in any browser. Controls: arrows/space, O overview, F fullscreen, S notes, ? help.`);

  if (args.open) {
    const opener = process.platform === 'win32' ? 'cmd' : process.platform === 'darwin' ? 'open' : 'xdg-open';
    const openArgs = process.platform === 'win32' ? ['/c', 'start', '', output] : [output];
    execFile(opener, openArgs, () => {});
  }
}

main().catch((e) => { console.error('✗ Conversion failed:', e.stack || e.message); process.exit(1); });
