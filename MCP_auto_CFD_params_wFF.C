R__LOAD_LIBRARY(/usr/opt/ddas/6.1-000/lib/libddaschannel.so)

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <TF1.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>

// --------------------------------------------
// Simple "does this trace fire?" helper
// - baseline from first Nbase samples
// - returns true if max deviation exceeds adcThr (ADC counts)
// --------------------------------------------
static bool TraceFired(const std::vector<uint16_t>& trace, double adcThr, int Nbase = 20)
{
  if (trace.empty()) return false;
  const int size = (int)trace.size();
  const int n = std::min(Nbase, size);

  double base = 0.0;
  for (int i = 0; i < n; i++) base += (double)trace[i];
  base /= (double)n;

  double maxDev = 0.0;
  for (int i = 0; i < size; i++) {
    double dev = std::abs((double)trace[i] - base);
    if (dev > maxDev) maxDev = dev;
  }
  return (maxDev > adcThr);
}

// ------------------------------------------------------------
// Pixie-like FAST FILTER (difference of two moving sums)
//
// ff[k] = sum_{i=0..L-1} y0[k-i]  -  sum_{i=0..L-1} y0[k-i-G]
//
// Parameters:
//   L = integration length (samples)
//   G = gap between the two sums (samples)
// ------------------------------------------------------------
static std::vector<double> FastFilter(const std::vector<double>& y0, int L, int G)
{
  const int n = (int)y0.size();
  std::vector<double> ff(n, 0.0);
  if (n <= 0) return ff;
  if (L <= 0) return ff;
  if (G <  0) G = 0;

  // prefix sum for O(1) window sums
  std::vector<double> ps(n + 1, 0.0);
  for (int i = 0; i < n; i++) ps[i + 1] = ps[i] + y0[i];

  auto sumRange = [&](int a, int b) -> double {
    // inclusive a..b
    if (a < 0) a = 0;
    if (b >= n) b = n - 1;
    if (a > b) return 0.0;
    return ps[b + 1] - ps[a];
  };

  // ff[k] uses windows ending at k
  for (int k = 0; k < n; k++) {
    int a1 = k - (L - 1);
    int b1 = k;
    int a2 = k - G - (L - 1);
    int b2 = k - G;

    // require both windows fully in range to avoid edge artifacts
    if (a2 < 0) continue;

    double s1 = sumRange(a1, b1);
    double s2 = sumRange(a2, b2);
    ff[k] = s1 - s2;
  }

  return ff;
}

// --------------------------------------------
// CFD t0 extractor: returns all t0s in a trace
// CFD is applied on FastFilter(y0)
// --------------------------------------------
static std::vector<double> GetT0sCFD(
    const std::vector<uint16_t>& trace,
    double f,
    int d,
    int NbaseEarly,
    int NbaseLate,
    int preGap,
    double negThr,   // used as +threshold magnitude (posThr = abs(negThr))
    int rearm,
    int FF_L,
    int FF_G
){
  std::vector<double> t0s;
  if (trace.empty()) return t0s;

  const int size = (int)trace.size();
  if (size <= d + 2) return t0s;

  // Convert to double
  std::vector<double> yd(size);
  for (int k = 0; k < size; k++) yd[k] = (double)trace[k];

  // Baseline from early window
  int nE = std::min(NbaseEarly, size);
  double baseEarly = 0.0;
  for (int k = 0; k < nE; k++) baseEarly += yd[k];
  baseEarly /= (double)nE;

  // Baseline-subtracted trace
  std::vector<double> y0(size);
  for (int k = 0; k < size; k++) y0[k] = yd[k] - baseEarly;

  // Fast filter output
  std::vector<double> ff = FastFilter(y0, FF_L, FF_G);

  // CFD on fast-filter output: y_cfd[k] = f*ff[k+d] - ff[k]
  std::vector<double> y_cfd(size, 0.0);
  for (int k = 0; k < size; k++) {
    double ff_shift = (k + d < size) ? ff[k + d] : 0.0;
    y_cfd[k] = f * ff_shift - ff[k];
  }

  // Find most positive point (anchor for "late" baseline)
  int kMax = d;
  double yMax = y_cfd[d];
  for (int k = d; k < size; k++) {
    if (y_cfd[k] > yMax) { yMax = y_cfd[k]; kMax = k; }
  }

  // If late, recompute baseline just before kMax and rebuild ff + CFD
  if (kMax > (NbaseLate + preGap)) {
    int k0 = kMax - preGap;
    int start = std::max(0, k0 - NbaseLate);

    double baseLate = 0.0;
    int cnt = 0;
    for (int k = start; k < k0; k++) { baseLate += yd[k]; cnt++; }
    if (cnt > 0) baseLate /= (double)cnt;

    for (int k = 0; k < size; k++) y0[k] = yd[k] - baseLate;

    ff = FastFilter(y0, FF_L, FF_G);

    for (int k = 0; k < size; k++) {
      double ff_shift = (k + d < size) ? ff[k + d] : 0.0;
      y_cfd[k] = f * ff_shift - ff[k];
    }
  }

  // Multi-pulse finder: require positive lobe then FALLING zero-crossing
  const double posThr = std::abs(negThr);
  bool sawPos = false;

  // Safe start index (avoid early FF zeros)
  int kStart = std::max(d + 1, FF_L + FF_G + d + 1);

  for (int k = kStart; k < size; k++) {
    if (y_cfd[k] > posThr) sawPos = true;

    if (sawPos && y_cfd[k-1] > 0.0 && y_cfd[k] <= 0.0) {
      double x1 = (double)(k - 1);
      double x2 = (double)k;
      double y1 = y_cfd[k-1];
      double y2 = y_cfd[k];

      double t0 = x1 + (0.0 - y1) * (x2 - x1) / (y2 - y1);
      t0s.push_back(t0);

      sawPos = false;
      k += rearm;
    }
  }

  return t0s;
}

struct BestResult {
  double sigma = 1e99;
  double f     = -1.0;
  int    d     = -1;
  double mean  = 0.0;
  Long64_t entries = 0;

  // new “quality” metrics
  double B = 0.0;            // flat background level from fit
  double purityFit = 0.0;    // fit-based purity in window
  double sideFrac  = 1.0;    // data-based fraction outside window
};

// ------------------------------
// Main macro
// ------------------------------
void MCP_auto_CFD_params_wFF()
{
  // ---------- Input ----------
  const char* infile = "/mnt/analysis/e20008/dumpedfiles/run2158-00.root";

  // Detector Upstream timing (crate/slot/channel)
  const int crateA = 1;
  const int slotA  = 2;
  const int chanA  = 1;

  // Detector Downstream timing (crate/slot/channel)
  const int crateB = 1;
  const int slotB  = 2;
  const int chanB  = 2;   // <-- EDIT if needed

  // ---------- OPTIONAL: Side gating on corner channels ----------
  const bool useSideGating = false;   // <-- set true to enable gating

  // Define the 4 corner channels for each detector (FILL THESE IN!)
  const int chanA_UL = 6;
  const int chanA_UR = 7;
  const int chanA_LL = 4;
  const int chanA_LR = 5;

  const int chanB_UL = 10;
  const int chanB_UR = 11;
  const int chanB_LL = 8;
  const int chanB_LR = 9;

  const bool requireSameSide = true;  // true: (A_L & B_L) OR (A_R & B_R)
  const double corner_adcThr = 500.0; // <-- tune

  // ---------- Baseline knobs ----------
  const int NbaseEarly = 25;
  const int NbaseLate  = 25;
  const int preGap     = 10;

  // ---------- CFD lobe threshold + rearm ----------
  const double negThr = 200.0; // used as +threshold magnitude (posThr = abs(negThr))
  const int rearm = 10;

  // ---------- Fast filter knobs (YOU TUNE THESE) ----------
  const int FF_L = /* YOUR_INPUT */ 8;
  const int FF_G = /* YOUR_INPUT */ 2;

  // If you want ns instead of samples
  const double dt_per_sample = 4.0;  // 250 MS/s -> 4 ns/sample

  // --------------------------
  // SCAN GRID (edit these)
  // --------------------------
  const double fMin = 0.40, fMax = 1.0, fStep = 0.025;
  const int    dMin = 1,    dMax = 6,   dStep = 1;

  // Histogram range (ns)
  const int nbins = 500;
  const double hmin = 50.0, hmax = 150.0;
  const double hmin_fit = 101.0, hmax_fit = 105.0;

  // ---------------------------------------------------------
  // New: quality cuts (tune these)
  // ---------------------------------------------------------
  const Long64_t minEntriesForFit = 50;   // reject low-stat histos
  const double sigmaMin = 0.10;           // reject pathological tiny sigma
  const double sigmaMax = 50.0;           // safety

  // Window around peak used to define “core” for purity metrics:
  // we use mean ± (coreNSigma * sigma)
  const double coreNSigma = 5.0;

  // Require cleanliness:
  const double purityMin = 0.70;          // fit-based signal fraction inside window
  const double sideFracMax = 0.20;        // data fraction outside window (scatter proxy)

  const int Nbest = 10; // keep top N

  // ---- Display / pause
  const int pause_ms = 250;

  // ---------- Open file / tree ----------
  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open file " << infile << "\n";
    return;
  }

  TTree* t = nullptr;
  f->GetObject("dchan", t);
  if (!t) {
    std::cerr << "ERROR: cannot find tree 'dchan'\n";
    return;
  }

  DDASEvent* ev = new DDASEvent();
  t->SetBranchAddress("ddasevent", &ev);

  // ---------------------------------------------------------
  // Cache traces AND coarse times once (so scan is fast)
  // ---------------------------------------------------------
  struct CachedEvent {
    std::vector<uint16_t> tA, tB;
    ULong64_t tsA = 0;
    ULong64_t tsB = 0;

    std::vector<uint16_t> A_UL, A_UR, A_LL, A_LR;
    std::vector<uint16_t> B_UL, B_UR, B_LL, B_LR;
  };

  std::vector<CachedEvent> events;
  events.reserve(t->GetEntries());

  const Long64_t nEntries = t->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";

  for (Long64_t i = 0; i < nEntries; i++) {
    t->GetEntry(i);

    CachedEvent ce;

    for (int j = 0; j < ev->GetNEvents(); j++) {
      ddaschannel* ch = ev->GetData()[j];
      if (!ch) continue;

      const int cr = ch->GetCrateID();
      const int sl = ch->GetSlotID();
      const int cn = ch->GetChannelID();

      if (cr == crateA && sl == slotA && cn == chanA) {
        ce.tA  = ch->GetTrace();
        ce.tsA = (ULong64_t)ch->GetCoarseTime();
      }

      if (cr == crateB && sl == slotB && cn == chanB) {
        ce.tB  = ch->GetTrace();
        ce.tsB = (ULong64_t)ch->GetCoarseTime();
      }

      if (useSideGating) {
        if (cn == chanA_UL && cr == crateA && sl == slotA) ce.A_UL = ch->GetTrace();
        if (cn == chanA_UR && cr == crateA && sl == slotA) ce.A_UR = ch->GetTrace();
        if (cn == chanA_LL && cr == crateA && sl == slotA) ce.A_LL = ch->GetTrace();
        if (cn == chanA_LR && cr == crateA && sl == slotA) ce.A_LR = ch->GetTrace();

        if (cn == chanB_UL && cr == crateB && sl == slotB) ce.B_UL = ch->GetTrace();
        if (cn == chanB_UR && cr == crateB && sl == slotB) ce.B_UR = ch->GetTrace();
        if (cn == chanB_LL && cr == crateB && sl == slotB) ce.B_LL = ch->GetTrace();
        if (cn == chanB_LR && cr == crateB && sl == slotB) ce.B_LR = ch->GetTrace();
      }
    }

    if (ce.tA.empty() || ce.tB.empty()) continue;

    if (useSideGating) {
      if (ce.A_UL.empty() || ce.A_UR.empty() || ce.A_LL.empty() || ce.A_LR.empty()) continue;
      if (ce.B_UL.empty() || ce.B_UR.empty() || ce.B_LL.empty() || ce.B_LR.empty()) continue;
    }

    events.emplace_back(std::move(ce));
  }

  std::cout << "Cached coincident events (timing present): " << events.size() << "\n";
  std::cout << "FF: L=" << FF_L << "  G=" << FF_G << "\n";
  std::cout << "Quality cuts: purityMin=" << purityMin
            << "  sideFracMax=" << sideFracMax
            << "  coreNSigma=" << coreNSigma << "\n";

  // --------------------------
  // Scan all (f,d) combinations
  // --------------------------
  std::cout << "\n=== Scan results ===\n";

  TCanvas* cScan = new TCanvas("cScan", "FF+CFD scan", 900, 600);

  // Fit function: gaussian + constant background
  //TF1 fSigBg("fSigBg", "gaus(0)+pol0(3)", hmin, hmax);
  //fSigBg.SetNpx(400);
  
  TF1 fGaus("fGaus", "gaus", hmin, hmax);
  fGaus.SetNpx(400);


  std::vector<BestResult> bestList;
  bestList.reserve(Nbest);

  // for annotation
  static TPaveText* pt = nullptr;

  for (double fTry = fMin; fTry <= fMax + 1e-12; fTry += fStep) {
    for (int dTry = dMin; dTry <= dMax; dTry += dStep) {

      TH1D hDt("hDt_tmp",
               Form("Delta t;#Delta t (ns);Counts (f=%.3f, d=%d, FF_L=%d, FF_G=%d)",
                    fTry, dTry, FF_L, FF_G),
               nbins, hmin, hmax);

      // Fill dt: tAbs = coarseTime + t0_ns
      for (const auto& ce : events) {

        if (useSideGating) {
          const bool A_left  = TraceFired(ce.A_UL, corner_adcThr) && TraceFired(ce.A_LL, corner_adcThr);
          const bool A_right = TraceFired(ce.A_UR, corner_adcThr) && TraceFired(ce.A_LR, corner_adcThr);

          const bool B_left  = TraceFired(ce.B_UL, corner_adcThr) && TraceFired(ce.B_LL, corner_adcThr);
          const bool B_right = TraceFired(ce.B_UR, corner_adcThr) && TraceFired(ce.B_LR, corner_adcThr);

          if (requireSameSide) {
            if (!((A_left && B_left) || (A_right && B_right))) continue;
          } else {
            if (!( (A_left || A_right) && (B_left || B_right) )) continue;
          }
        }

        auto t0A = GetT0sCFD(ce.tA, fTry, dTry, NbaseEarly, NbaseLate, preGap, negThr, rearm, FF_L, FF_G);
        auto t0B = GetT0sCFD(ce.tB, fTry, dTry, NbaseEarly, NbaseLate, preGap, negThr, rearm, FF_L, FF_G);
        if (t0A.empty() || t0B.empty()) continue;

        const double t0A_ns = t0A.front() * dt_per_sample;
        const double t0B_ns = t0B.front() * dt_per_sample;

        const double tAbsA = (double)ce.tsA + t0A_ns;
        const double tAbsB = (double)ce.tsB + t0B_ns;

        const double dt = tAbsB - tAbsA;
        hDt.Fill(dt);
      }

      const Long64_t Ntot = (Long64_t)hDt.GetEntries();
      if (Ntot < minEntriesForFit) {
        std::cout << "f=" << fTry << " d=" << dTry
                  << "  entries=" << Ntot
                  << "  SKIP (low stats)\n";
        continue;
      }

      // peak guess
      int maxBin = hDt.GetMaximumBin();
      double xPeak = hDt.GetXaxis()->GetBinCenter(maxBin);
      double yPeak = hDt.GetBinContent(maxBin);

      // Fit over full range but with good initial params
      fGaus.SetRange(hmin_fit, hmax_fit);
      fGaus.SetParameters(std::max(1.0, yPeak), xPeak, 2.0);//, 0.0); // A, mean, sigma, B
      hDt.Fit(&fGaus, "R"); // quiet, no draw

      double A     = fGaus.GetParameter(0);
      double mean  = fGaus.GetParameter(1);
      double sigma = std::abs(fGaus.GetParameter(2));
      //double B     = fSigBg.GetParameter(3);
      //if (B < 0) B = 0.0;

      // sanity
      if (!(sigma > sigmaMin && sigma < sigmaMax)) {
        std::cout << "f=" << fTry << " d=" << dTry
                  << "  entries=" << Ntot
                  << "  sigma=" << sigma
                  << "  REJECT (sigma out of bounds)\n";
        continue;
      }

      // ---- Define “core” window around fitted mean ----
      double win = coreNSigma * sigma;
      // keep it inside hist range
      double x1 = std::max(hmin, mean - win);
      double x2 = std::min(hmax, mean + win);
      
      // --- Data-based background estimate from sidebands (no background fit) ---
      // Sidebands = everything outside [x1,x2]
      int bL1 = 1;
      int bL2 = hDt.GetXaxis()->FindBin(x1) - 1;
      int bR1 = hDt.GetXaxis()->FindBin(x2) + 1;
      int bR2 = hDt.GetNbinsX();

      double Nside = 0.0;
      int NsideBins = 0;

      if (bL2 >= bL1) { Nside += hDt.Integral(bL1, bL2); NsideBins += (bL2 - bL1 + 1); }
      if (bR2 >= bR1) { Nside += hDt.Integral(bR1, bR2); NsideBins += (bR2 - bR1 + 1); }

      double B = (NsideBins > 0) ? (Nside / (double)NsideBins) : 0.0; // counts/bin

      // ---- Fit-based purity estimate: gaussian area vs background area in [x1,x2] ----
      const double binw = hDt.GetXaxis()->GetBinWidth(1);

      // Gaussian area in [x1,x2] -> counts
      TF1 gtmp("gtmp","gaus", x1, x2);
      gtmp.SetParameters(A, mean, sigma);
      const double S  = gtmp.Integral(x1, x2) / binw;

      // Background counts in [x1,x2] from sideband estimate
      const double Bg = B * ((x2 - x1) / binw);

      const double purityFit = (S > 0.0) ? (S / (S + Bg)) : 0.0;

      // Data-based scatter proxy remains the same:
      int b1 = hDt.GetXaxis()->FindBin(x1);
      int b2 = hDt.GetXaxis()->FindBin(x2);
      double Nwin = hDt.Integral(b1, b2);
      double sideFrac = (Ntot > 0) ? (1.0 - Nwin / (double)Ntot) : 1.0;
      
      // ---- Apply cleanliness cuts ----
      const bool passQuality = (purityFit >= purityMin) && (sideFrac <= sideFracMax);

      // ---- Draw + overlay the fit for the user ----
      cScan->cd();
      hDt.SetLineColor(kBlue+1);
      hDt.Draw();

      // draw fit on top (now we DO draw)
      hDt.Fit(&fGaus, "R");   // not quiet, draws the curve
      fGaus.SetLineColor(kRed);

      if (pt) { delete pt; pt = nullptr; }
      pt = new TPaveText(0.12, 0.65, 0.42, 0.95, "NDC");
      pt->SetFillColor(0);
      pt->SetBorderSize(1);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.028);

      pt->AddText(Form("f = %.3f, d = %d", fTry, dTry));
      pt->AddText(Form("entries = %lld", Ntot));
      pt->AddText(Form("#mu = %.3f ns", mean));
      pt->AddText(Form("#sigma = %.3f ns  (FWHM=%.3f ns)", sigma, 2.355*sigma));
      pt->AddText(Form("B = %.2f (cts/bin)", B));
      pt->AddText(Form("core = #pm %.1f#sigma  -> [%.2f, %.2f] ns", coreNSigma, x1, x2));
      pt->AddText(Form("purityFit = %.3f", purityFit));
      pt->AddText(Form("sideFrac(data) = %.3f", sideFrac));
      pt->AddText(passQuality ? "QUALITY: PASS" : "QUALITY: FAIL");
      pt->Draw("SAME");

      cScan->Modified();
      cScan->Update();
      gSystem->ProcessEvents();
      gSystem->Sleep(pause_ms);

      // ---- Print one-line summary ----
      std::cout << "f=" << fTry << " d=" << dTry
                << "  entries=" << Ntot
                << "  mean=" << mean
                << "  sigma=" << sigma
                << "  purityFit=" << purityFit
                << "  sideFrac=" << sideFrac
                << (passQuality ? "  PASS\n" : "  FAIL\n");

      if (!passQuality) continue;

      // ---- Keep top N best among PASS solutions (rank by sigma) ----
      BestResult cand;
      cand.sigma = sigma;
      cand.f = fTry;
      cand.d = dTry;
      cand.mean = mean;
      cand.entries = Ntot;
      cand.B = B;
      cand.purityFit = purityFit;
      cand.sideFrac = sideFrac;

      if ((int)bestList.size() < Nbest || cand.sigma < bestList.back().sigma) {
        bestList.push_back(cand);
        std::sort(bestList.begin(), bestList.end(),
                  [](const BestResult& a, const BestResult& b) {
                    return a.sigma < b.sigma;
                  });
        if ((int)bestList.size() > Nbest) bestList.pop_back();
      }
    }
  }

  std::cout << "\n=== TOP " << bestList.size() << " RESULTS (PASS quality, smallest sigma) ===\n";
  for (size_t i = 0; i < bestList.size(); i++) {
    const auto& r = bestList[i];
    std::cout << i+1 << ") "
              << "sigma=" << r.sigma << " ns"
              << "  FWHM=" << 2.355 * r.sigma << " ns"
              << "  purityFit=" << r.purityFit
              << "  sideFrac=" << r.sideFrac
              << "  B=" << r.B
              << "  f=" << r.f
              << "  d=" << r.d
              << "  mean=" << r.mean
              << "  entries=" << r.entries
              << "\n";
  }

  std::cout << "\nDone.\n";
}
