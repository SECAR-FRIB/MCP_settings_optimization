R__LOAD_LIBRARY(/usr/opt/ddas/6.1-000/lib/libddaschannel.so)
//runs with root6.24/06

#include <TSystem.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <iostream>
#include <sstream>
#include <vector>


// ------------------------------
// Simple trace-plotting function
// ------------------------------
void MCP_timing_trace_analysis() {
    
    // --------------------------
    // User knobs
    // --------------------------
    const double f =0.925;     // scale factor for inverted signal
    const int d  = 4;     // delay in samples (can be non-integer too)

    // Open ROOT file
    TFile* pFile = new TFile("/pathTo/dumpedfiles/run2158-00.root");
    // Choose: Crate / Module / Channel
    int plotCrate = 1; // Use 0 for target and 1 for FP4
    int plotMod = 0+2; // Module id's 0 and 1 are reserved for the computer in the root scripts
    int plotChan = 2;
    // How many traces to draw
    const int maxTracesToPlot =1;
    
    if (!pFile || pFile->IsZombie()) 
    {
      std::cerr << "ERROR: Could not open file." << std::endl;
      return;
    }
    
    // Get tree
    TTree* pTree = nullptr;
    pFile->GetObject("dchan", pTree);
    if (!pTree) 
    {
      std::cerr << "ERROR: Tree 'dchan' not found." << std::endl;
      return;
    }

    // Set up event object
    DDASEvent* pEvent = new DDASEvent();
    pTree->SetBranchAddress("ddasevent", &pEvent);
    
    //TCanvas* c1 = new TCanvas("c1", "Traces", 1500, 1000);
    TCanvas* c2 = new TCanvas("c2", "CFD", 1500, 1000);
   
    
    int tracesPlotted = 0;
    int cfdPlotted    = 0;

    // Predefined list of colors to cycle through
    int colors[] = {kBlue, kRed,  kGreen+2, kMagenta, kOrange+7, kCyan+2 };

    // Loop over entries in the TTree
    for (int i = 0; i < pTree->GetEntries(); ++i) 
    {  
      pTree->GetEntry(i);

      std::vector<uint16_t> trace;

      // Loop over channels in this event
      for (int j = 0; j < pEvent->GetNEvents(); j++) 
      {

        ddaschannel* dchan = pEvent->GetData()[j];

        // Select your desired crate/slot/channel
        if (dchan->GetCrateID() == plotCrate &&
            dchan->GetSlotID()  == plotMod &&
            dchan->GetChannelID() == plotChan)
        {
          trace = dchan->GetTrace();
        }
      }

      // Skip if no trace in this event
      if (trace.empty()) continue;

      int size = (int)trace.size();
      
      // Convert trace into doubles first
      std::vector<double> xd(size), yd(size);
      for (int k = 0; k < size; k++) {
        xd[k] = k;
        yd[k] = (double)trace[k];
      }
      
      const int NbaseEarly = 5;   // fallback if too early
      const int NbaseLate  = 25;   // baseline length for late pulses
      const int preGap     = 10;   // how far before the pulse to stop baseline window
      
      // baseline-subtract using early baseline FIRST (rough)
      double baseEarly = 0;
      for (int k=0; k<std::min(NbaseEarly,size); k++) baseEarly += yd[k];
      baseEarly /= (double)std::min(NbaseEarly,size);
      
      // Baseline-subtracted trace
      std::vector<double> y0(size);
      for (int k=0; k<size; k++) y0[k] = yd[k] - baseEarly;
            
      // CFD waveform: y(k-d) - f*y(k)
      std::vector<double> y_cfd(size, 0.0);
      for (int k = 0; k < size; k++) {
        double y_shift = (k + d < size) ? y0[k + d] : 0.0;
        y_cfd[k] = f * y_shift - y0[k];
      }
      
      // find positive peak index kMin (pulse anchor)
      int kMax = d;
      double yMax = y_cfd[d];
      for (int k = d; k < size; k++) {
        if (y_cfd[k] > yMax) { yMax = y_cfd[k]; kMax = k; }
      }
      
      // if the pulse is "late", recompute baseline using a window before kMax
      if (kMax > (NbaseLate + preGap)) {
        double baseLate = 0.0;
        int k0 = kMax - preGap;                // end of baseline window
        int start = std::max(0, k0 - NbaseLate);
        for (int k = start; k < k0; k++) baseLate += yd[k];
        baseLate /= (double)NbaseLate;
      
        // recompute y0 and CFD with the better baseline
        for (int k=0; k<size; k++) y0[k] = yd[k] - baseLate;
      
        for (int k=0; k<size; k++) {
          double y_shift = (k + d < size) ? y0[k + d] : 0.0;
          y_cfd[k] = f * y_shift - y0[k];
        }
      }


      // ---------- Draw original trace ----------
      TGraph* gr = new TGraph(size, xd.data(), yd.data());
      int color = colors[tracesPlotted % (sizeof(colors)/sizeof(int))];
      gr->SetLineColor(color);
      gr->SetLineWidth(1);
      gr->SetLineStyle(2);
      
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(1);

      // Draw first trace normally, the rest "same"
      c2->cd();
      if (tracesPlotted == 0)
	    {
        // Build the title string
        std::stringstream ss;
        ss <<  maxTracesToPlot << " Traces of"
           << " Mod " << plotMod << " Chan " << plotChan;
        gr->SetTitle(ss.str().c_str());
        gr->Draw("ALP");
        //gr->GetXaxis()->SetRangeUser(-100, 100);
        gr->GetYaxis()->SetRangeUser(-25000, 25000);
      }
      else
        gr->Draw("L SAME");

      tracesPlotted++;
      
      
      // ---------- Draw CFD ----------
      TGraph* grCFD = new TGraph(size, xd.data(), y_cfd.data());
      grCFD->SetLineColor(color);
      grCFD->SetLineWidth(2);

      // Title + draw
      c2->cd();
      if (cfdPlotted == 0)
      {
        grCFD->SetTitle(Form("CFD: y(k-%d) - %.3f*y(k)", d, f));
        grCFD->Draw("L SAME");
        
        // dashed y=0 line across the current pad x-range
        gPad->Update();
        double xmin = gPad->GetUxmin();
        double xmax = gPad->GetUxmax();
        TLine* l0 = new TLine(xmin, 0.0, xmax, 0.0);
        l0->SetLineStyle(2);   // dashed
        l0->SetLineWidth(2);
        l0->SetLineColor(kGray+2);
        l0->Draw("SAME");
      } else {
        grCFD->Draw("L SAME");
      }
      cfdPlotted++;
      
      // ============================================================
      // Multi-pulse t0 finder (can return 0, 1, 2... t0 values)
      // Sequence: must go below thresh, then take next rising crossing
      // ============================================================
      
      std::vector<double> t0s;

      const double thresh = 400.0;   // tune: require CFD to go below this to "arm"
      const int rearm = 10;           // tune: skip this many samples after a found crossing
      const double settleThr = 200.0; // tune: must settle near 0 before re-arming (optional)

      bool armed = true;

      bool sawPos = false;
      bool sawNeg = false;
      
      for (int k = d + 1; k < size; k++)
      {
        if (armed)
        {
          //if (y_cfd[k] < thresh) sawNeg = true;
          if (y_cfd[k] > thresh) sawPos = true;
          //if (sawNeg && y_cfd[k-1] < 0.0 && y_cfd[k] >= 0.0)
          if (sawPos && y_cfd[k-1] > 0.0 && y_cfd[k] <= 0.0)
          {
            // linear interpolation for t0
            double x1 = xd[k-1], x2 = xd[k];
            double y1 = y_cfd[k-1], y2 = y_cfd[k];
            double t0 = x1 + (0.0 - y1) * (x2 - x1) / (y2 - y1);

            t0s.push_back(t0);

            // mark on plot
            TMarker* m0 = new TMarker(t0, 0.0, 20);
            m0->SetMarkerColor(kRed);
            m0->SetMarkerSize(1.2);
            m0->Draw("SAME");

            // disarm briefly
            armed = false;
            sawNeg = false;
            sawPos = false;
            k += rearm;
          }
        }
        else
        {
          // rearm once we're close to baseline again (helps avoid ringing crossings)
          armed = true;
        }
      }

      std::cout << "Event " << i << " found " << t0s.size() << " t0(s): ";
      for (double v : t0s) std::cout << v << " ";
      std::cout << std::endl;
    
      for (size_t i0 = 0; i0 < t0s.size(); i0++) {
        std::cout << "kMax=" << kMax
                  << " yMax=" << yMax
                  << " t0[" << i0 << "]=" << t0s[i0]
                  << std::endl;
      }



      if (tracesPlotted >= maxTracesToPlot)
        break;
    }
    
    std::cout << "Plotted " << tracesPlotted << " traces." << std::endl;
}
