R__LOAD_LIBRARY(/usr/opt/ddas/6.1-000/lib/libddaschannel.so)

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>

// ------------------------------------------------------------
// Pixie-like FAST FILTER (difference of two moving sums)
//
// ff[k] = sum_{i=0..L-1} y0[k-i]  -  sum_{i=0..L-1} y0[k-i-G]
//
// Parameters:
//   L = integration length (samples)
//   G = gap between the two sums (samples)
//
// Notes:
// - returns vector same size as y0
// - valid output begins around k >= (L+G) ; earlier points are 0
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
// NOW: CFD is applied on FastFilter(y0)
// --------------------------------------------
static std::vector<double> GetT0sCFD(
    const std::vector<uint16_t>& trace,
    double f,
    int d,
    int NbaseEarly,
    int NbaseLate,
    int preGap,
    double negThr,
    int rearm,
    // --- Fast filter knobs (Pixie-like) ---
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

  // --- Fast filter output ---
  // ff is what Pixie uses for timing/CFD-like operations (conceptually)
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

  // Multi-pulse finder:
  // require a positive lobe, then take the FALLING zero-crossing (+ -> -)
  const double posThr = std::abs(negThr);  // keep your knob semantics
  bool sawPos = false;

  // NOTE: fast filter has a natural "dead" region at the beginning.
  // A safe start index is: FF_L + FF_G + d + 1
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


// ------------------------------
// Main macro
// ------------------------------
void MCP_plot_dt_histo_ff_cfd()
{
  // ---------- Input ----------
  const char* infile = "/pathTO/dumpedfiles/run2158-00.root";

  // Detector A (crate/slot/channel)
  const int crateA = 1;
  const int slotA  = 2;
  const int chanA  = 1;   // Upstream MCP

  // Detector B (crate/slot/channel)
  const int crateB = 1;
  const int slotB  = 2;
  const int chanB  = 2;   // Downstream MCP

  // ---------- CFD knobs ----------
  const double fA = 0.925;
  const int    dA = 4;

  const double fB = fA;
  const int    dB = dA;

  const int NbaseEarly = 25;
  const int NbaseLate  = 25;
  const int preGap     = 10;

  const double negThr = 200.0; // (used as +threshold magnitude here)
  const int rearm = 10;

  // ---------- Fast filter knobs (YOU TUNE THESE) ----------
  // Typical starting points (250 MS/s): L=4..8, G=1..4
  const int FF_L = /* YOUR_INPUT (e.g. 4) */ 8;
  const int FF_G = /* YOUR_INPUT (e.g. 2) */ 2;

  // If you want ns instead of samples
  const double dt_per_sample = 4.0;  // 250 MS/s -> 4 ns

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

  // ---------- Histogram ----------
  TH1D* hDt = new TH1D("hDt","Delta t = t0(DMCP) - t0(UMCP);#Deltat (0.2 ns bin); Counts",
                      500, 50, 150);

  // ---------- Loop events ----------
  const Long64_t nEntries = t->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";

  for (Long64_t i = 0; i < nEntries; i++) {
    t->GetEntry(i);

    std::vector<uint16_t> traceA, traceB;
    ULong64_t tsA = 0, tsB = 0;
    ULong64_t tA  = 0, tB  = 0;
    bool foundA = false, foundB = false;

    for (int j = 0; j < ev->GetNEvents(); j++) {
      ddaschannel* ch = ev->GetData()[j];
      if (!ch) continue;

      const int cr = ch->GetCrateID();
      const int sl = ch->GetSlotID();
      const int cn = ch->GetChannelID();

      if (cr == crateA && sl == slotA && cn == chanA) {
        traceA = ch->GetTrace();
        tsA    = (ULong64_t)ch->GetCoarseTime(); // in ns (per DDAS)
        tA     = (ULong64_t)ch->GetTime();       // in ns (HW timing)
        foundA = true;
      }

      if (cr == crateB && sl == slotB && cn == chanB) {
        traceB = ch->GetTrace();
        tsB    = (ULong64_t)ch->GetCoarseTime(); // in ns
        tB     = (ULong64_t)ch->GetTime();       // in ns
        foundB = true;
      }
    }

    if (!foundA || !foundB) continue;
    if (traceA.empty() || traceB.empty()) continue;

    auto t0A = GetT0sCFD(traceA, fA, dA, NbaseEarly, NbaseLate, preGap, negThr, rearm, FF_L, FF_G);
    auto t0B = GetT0sCFD(traceB, fB, dB, NbaseEarly, NbaseLate, preGap, negThr, rearm, FF_L, FF_G);
    if (t0A.empty() || t0B.empty()) continue;

    // Convert t0 (samples) to ns
    // NOTE: remove/adjust the "+8" once you understand your kRef/alignment
    const double t0A_ns = t0A.front() * dt_per_sample;
    const double t0B_ns = t0B.front() * dt_per_sample;

    // Coarse time is already in ns (DDAS GetCoarseTime)
    const double tsA_ns = (double)tsA;
    const double tsB_ns = (double)tsB;

    // Your convention (keep as you had it)
    const double tAbsA = tsA_ns + t0A_ns;
    const double tAbsB = tsB_ns + t0B_ns;

    const double dt = tAbsB - tAbsA; // ns
    hDt->Fill(dt);
    //if(dt<100) {std::cout<<i<<"\t"<<t0A.front()<<"\t"<<t0B.front()<<"\t"<<tB-tA<<"\t"<<dt<<endl;}
    
    // Hardware dt (for comparison)
    //const double dt_hw = (double)tB - (double)tA;
    //hDt->Fill(dt_hw);
  }

  // ---------- Save + Draw ----------
  TFile* out = TFile::Open("MCP_dt_hist.root", "RECREATE");
  gStyle->SetOptStat(1000111);
  hDt->Write();
  out->Close();

  TCanvas* c = new TCanvas("c_dt", "Delta t", 900, 600);
  hDt->Draw();
  c->Update();
  std::cout << "bin width " << hDt->GetXaxis()->GetBinWidth(1) << " ns\n";

  std::cout << "Wrote histogram to MCP_dt_hist.root\n";
}
