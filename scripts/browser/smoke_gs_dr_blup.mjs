/**
 * Puppeteer smoke: Growth Summaries → Trait DR → BLUP (+ Clippy alive).
 *
 * Prerequisites: StatFaRmer main app already listening, soy project loaded.
 *
 *   node scripts/browser/smoke_gs_dr_blup.mjs --url http://127.0.0.1:3848
 *
 * Artifacts: outputs/measurements/puppeteer/smoke_*
 */
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import puppeteer from "puppeteer";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(__dirname, "../..");
const outDir = path.join(root, "outputs", "measurements", "puppeteer");

function arg(name, fallback) {
  const i = process.argv.indexOf(name);
  if (i >= 0 && process.argv[i + 1]) return process.argv[i + 1];
  return fallback;
}

const url = arg("--url", process.env.STATFARMER_URL || "http://127.0.0.1:3848");
const settleMs = Number(arg("--settle", "2500"));

fs.mkdirSync(outDir, { recursive: true });

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function clickTabByText(page, re) {
  return page.evaluate((pattern) => {
    const rx = new RegExp(pattern, "i");
    const links = Array.from(
      document.querySelectorAll(".nav-tabs a, .nav-tabs .nav-link, a[data-toggle='tab'], a[data-bs-toggle='tab']")
    );
    const hit = links.find((a) => rx.test((a.textContent || "").trim()));
    if (!hit) return { ok: false, texts: links.map((a) => (a.textContent || "").trim()).slice(0, 30) };
    hit.click();
    return { ok: true, text: (hit.textContent || "").trim() };
  }, re.source || re);
}

async function clickButtonById(page, id) {
  // Shiny ns uses hyphen: growth_summaries-calculate_growth
  const safe = `#${id.replace(/([:.])/g, "\\$1")}`;
  const el = await page.$(safe);
  if (!el) {
    // Fallback: click by evaluating id (avoids CSS.escape quirks)
    return page.evaluate((buttonId) => {
      const el = document.getElementById(buttonId);
      if (!el) return false;
      el.click();
      return true;
    }, id);
  }
  await el.click();
  return true;
}

async function shinyBusyGone(page, timeoutMs = 120000) {
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    const busy = await page.evaluate(() => {
      const b = document.getElementById("shiny-busy") || document.querySelector(".shiny-busy");
      if (b) return true;
      const html = document.documentElement;
      return html && html.classList.contains("shiny-busy");
    });
    if (!busy) return true;
    await sleep(400);
  }
  return false;
}

const report = {
  url,
  started_at: new Date().toISOString(),
  steps: [],
  console_errors: [],
  page_errors: [],
  ok: false,
};

function resolveChrome() {
  if (process.env.PUPPETEER_EXECUTABLE_PATH) return process.env.PUPPETEER_EXECUTABLE_PATH;
  const candidates = [
    "C:\\\\Program Files\\\\Google\\\\Chrome\\\\Application\\\\chrome.exe",
    "C:\\\\Program Files (x86)\\\\Google\\\\Chrome\\\\Application\\\\chrome.exe",
    "C:\\\\Program Files\\\\Microsoft\\\\Edge\\\\Application\\\\msedge.exe",
    "C:\\\\Program Files (x86)\\\\Microsoft\\\\Edge\\\\Application\\\\msedge.exe",
  ];
  for (const p of candidates) {
    if (fs.existsSync(p)) return p;
  }
  return undefined;
}

const chromePath = resolveChrome();
const browser = await puppeteer.launch({
  headless: true,
  executablePath: chromePath,
  args: ["--no-sandbox", "--disable-setuid-sandbox"],
});
report.chrome = chromePath || "(puppeteer default)";
const page = await browser.newPage();
page.setDefaultTimeout(90000);
await page.setViewport({ width: 1500, height: 960 });

page.on("pageerror", (err) => {
  report.page_errors.push(String(err && err.message ? err.message : err));
});
page.on("console", (msg) => {
  if (msg.type() === "error") {
    report.console_errors.push(msg.text());
  }
});

function pushStep(name, extra = {}) {
  report.steps.push({ name, at: new Date().toISOString(), ...extra });
}

try {
  await page.goto(url, { waitUntil: "networkidle2" });
  await page.waitForSelector("body");
  await sleep(settleMs);
  await shinyBusyGone(page, 60000);
  await page.screenshot({ path: path.join(outDir, "smoke_01_landing.png"), fullPage: true });
  pushStep("landing");

  // --- Growth Summaries ---
  let tab = await clickTabByText(page, "Growth Summaries");
  pushStep("open_growth_summaries", tab);
  if (!tab.ok) throw new Error(`Growth Summaries tab not found; saw: ${(tab.texts || []).join(" | ")}`);
  await sleep(settleMs);
  await shinyBusyGone(page);

  const gsClicked = await clickButtonById(page, "growth_summaries-calculate_growth");
  pushStep("click_calculate_summaries", { clicked: gsClicked });
  if (!gsClicked) throw new Error("Calculate Summaries button not found");
  await sleep(2000);
  await shinyBusyGone(page, 180000);
  await page.screenshot({ path: path.join(outDir, "smoke_02_growth.png"), fullPage: true });

  // Wait for DT or at least no shiny notification of hard error
  const gsOk = await page.evaluate(() => {
    const err = document.querySelector(".shiny-output-error");
    const table = document.querySelector("#growth_summaries-growth_summary_table");
    const note = document.querySelector("#growth_summaries-summary_note");
    const body = (document.body.innerText || "");
    return {
      hasError: !!(err && err.offsetParent !== null && /error|validate|unused/i.test(err.textContent || "")),
      errorText: err ? (err.textContent || "").slice(0, 300) : "",
      hasTable: !!table,
      hasNamedMetrics:
        /__auc|named metrics|trait__/.test(body) ||
        /__auc|named metrics/.test((note && note.textContent) || ""),
    };
  });
  pushStep("growth_result", gsOk);
  if (gsOk.hasError) throw new Error(`Growth Summaries UI error: ${gsOk.errorText}`);
  if (!gsOk.hasNamedMetrics) {
    throw new Error("Growth Summaries did not expose trait__metric / named metrics text");
  }

  // --- Trait DR ---
  tab = await clickTabByText(page, "Trait DR|Trait Reduction|Dimensionality");
  pushStep("open_trait_dr", tab);
  if (!tab.ok) {
    // try exact known title
    tab = await clickTabByText(page, "Trait");
  }
  if (!tab.ok) throw new Error(`Trait DR tab not found; saw: ${(tab.texts || []).join(" | ")}`);
  await sleep(settleMs);
  await shinyBusyGone(page);

  // Select first two traits in selectize if present
  await page.evaluate(() => {
    const sel = document.querySelector("#trait_reduction-dr_traits");
    if (!sel) return;
    const opts = Array.from(sel.options || []).map((o) => o.value).filter(Boolean);
    if (opts.length >= 2 && window.Shiny) {
      Shiny.setInputValue("trait_reduction-dr_traits", opts.slice(0, Math.min(6, opts.length)), {
        priority: "event",
      });
    }
  });
  await sleep(800);

  const drClicked = await clickButtonById(page, "trait_reduction-calculate_dr");
  pushStep("click_calculate_dr", { clicked: drClicked });
  if (drClicked) {
    await sleep(2000);
    await shinyBusyGone(page, 180000);
  }
  await page.screenshot({ path: path.join(outDir, "smoke_03_trait_dr.png"), fullPage: true });

  const drErr = await page.evaluate(() => {
    const nodes = Array.from(document.querySelectorAll(".shiny-output-error, .shiny-output-error-validation"));
    return nodes
      .map((n) => (n.textContent || "").trim())
      .filter((t) => t.length)
      .slice(0, 5);
  });
  pushStep("trait_dr_errors", { errors: drErr });
  if (drErr.some((t) => /unused argument|неиспользован|jsonlite/i.test(t))) {
    throw new Error(`Trait DR validate mask still broken: ${drErr.join(" || ")}`);
  }

  // --- BLUP ---
  tab = await clickTabByText(page, "^BLUP|BLUP Ranking");
  pushStep("open_blup", tab);
  if (!tab.ok) tab = await clickTabByText(page, "BLUP");
  if (!tab.ok) throw new Error(`BLUP tab not found; saw: ${(tab.texts || []).join(" | ")}`);
  await sleep(settleMs);
  await shinyBusyGone(page);

  const blupClicked = await clickButtonById(page, "blup_ranking-calculate_blup");
  pushStep("click_calculate_blup", { clicked: blupClicked });
  if (!blupClicked) throw new Error("Calculate BLUP button not found");
  await sleep(2000);
  await shinyBusyGone(page, 180000);
  await page.screenshot({ path: path.join(outDir, "smoke_04_blup.png"), fullPage: true });

  const blupProbe = await page.evaluate(() => {
    const nodes = Array.from(document.querySelectorAll(".shiny-output-error, .shiny-output-error-validation"));
    const errors = nodes.map((n) => (n.textContent || "").trim()).filter(Boolean).slice(0, 8);
    const clippy = document.querySelector("#farma_clip-clippy_root, #farma_clip-panel, .farma-clip, [id*='farma_clip']");
    const table = document.querySelector("#blup_ranking-blup_table, #blup_ranking-genotype_table");
    return {
      errors,
      hasClippy: !!clippy,
      hasTableish: !!table || !!document.querySelector("#blup_ranking-rank_table"),
      bodySnippet: (document.body.innerText || "").slice(0, 400),
    };
  });
  pushStep("blup_result", blupProbe);

  const bad = [...blupProbe.errors, ...report.page_errors, ...report.console_errors].filter((t) =>
    /unused argument|неиспользован|jsonlite::validate|Error in validate/i.test(t)
  );
  if (bad.length) {
    throw new Error(`BLUP/Clippy still hit validate mask: ${bad.join(" || ")}`);
  }

  // Clippy should not be a hard shiny error panel
  const clippyBroken = blupProbe.errors.some((t) => /farma_clip|clippy/i.test(t));
  if (clippyBroken) throw new Error("Clippy panel still erroring");

  report.ok = true;
} catch (err) {
  report.ok = false;
  report.error = String(err && err.message ? err.message : err);
  try {
    await page.screenshot({ path: path.join(outDir, "smoke_error.png"), fullPage: true });
  } catch (_) {
    /* ignore */
  }
} finally {
  report.finished_at = new Date().toISOString();
  fs.writeFileSync(path.join(outDir, "smoke_gs_dr_blup.json"), JSON.stringify(report, null, 2));
  await browser.close();
}

console.log(report.ok ? `OK → ${path.join(outDir, "smoke_gs_dr_blup.json")}` : `FAIL → ${report.error}`);
process.exit(report.ok ? 0 : 1);
