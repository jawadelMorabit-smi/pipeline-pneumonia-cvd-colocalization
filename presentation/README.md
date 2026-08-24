# pptx → interactive HTML

A small system that converts **any** PowerPoint `.pptx` into a single, self-contained
**interactive HTML presentation** — no server, no external assets, works offline, and is
trivially shareable (just send the one `.html` file).

## Usage

```bash
npm install                      # one-time: installs jszip + fast-xml-parser
node generate.js <input.pptx> [output.html] [options]
```

Examples:

```bash
node generate.js Pneumonia_CVD_Colocalization_deck.pptx
node generate.js deck.pptx talk.html --open
node generate.js deck.pptx --fragments --title "My Talk"
```

Default output is `<input>.interactive.html` next to the source file.

### Options

| Option        | Effect                                                            |
|---------------|-------------------------------------------------------------------|
| `--fragments` | Reveal bullets/shapes step-by-step on ← / → (default: whole slide)|
| `--title "…"` | Document title (default: derived from the file name)              |
| `--open`      | Open the result in your default browser when done                 |

## Controls (in the generated deck)

| Key / gesture                         | Action                    |
|---------------------------------------|---------------------------|
| → · Space · PageDown · click          | Next step / slide         |
| ← · PageUp                            | Previous                  |
| ↑ / ↓                                 | Previous / next slide     |
| Home · End                            | First / last slide        |
| number then Enter                     | Jump to slide             |
| **O**                                 | Overview grid (click a thumbnail to jump) |
| **F**                                 | Fullscreen                |
| **S**                                 | Speaker notes panel       |
| **?**                                 | Keyboard help             |
| Esc                                   | Close overlay / exit overview |
| Swipe                                 | Navigate (touch devices)  |
| Click an image                        | Zoom in / out             |

The current slide is reflected in the URL (`…#/12`), so any slide is bookmarkable and
shareable as a deep link. The stage auto-scales to any screen size.

## How it works

```
generate.js            CLI: args → parse → build → write
lib/pptx-parser.js     unzips the .pptx and reads the OOXML into a plain slide model:
                       slide order, geometry (EMU→px), text runs + formatting,
                       theme/​direct colours, images (inlined as data URIs), tables,
                       grouped shapes, and speaker notes.
lib/html-builder.js    renders the model to one HTML file — each slide is a fixed
                       1280×720 stage of absolutely-positioned elements — and injects
                       the interactive runtime (navigation, overview, notes, zoom…).
```

Because every element is absolutely positioned, fidelity to the original layout is high
and step-by-step reveals never reflow the slide (they only change opacity).

### Fidelity notes / limits
- Text, colours (direct + theme), fonts, sizes, bold/italic/underline, alignment,
  bullets with hanging indents, shape fills, rounded rectangles, borders, images and
  tables are reproduced.
- Fonts render with the viewer's locally-installed families (e.g. Calibri, Cambria);
  exotic fonts fall back gracefully.
- Charts/SmartArt saved as vector-only objects (no fallback image) and slide
  transitions/animations authored in PowerPoint are not reproduced — the deck gets its
  own clean transitions instead.
