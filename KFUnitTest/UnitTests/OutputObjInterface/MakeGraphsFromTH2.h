#pragma once

template <typename T, typename U> U *MakeIntegatedHistoFromTH2(T fHisto2D) {
  U *fHistOut = new U();
  fHistOut = fHisto2D->ProjectionY();
  ;
  return fHistOut;
}

template <typename T> TGraphErrors *MakeMeanGraphFromTH2(T fHisto2D) {
  int num_of_bins = fHisto2D->GetXaxis()->GetLast();
  float x[num_of_bins], y[num_of_bins], xerr[num_of_bins], yerr[num_of_bins];
  for (int ibin = 1; ibin <= num_of_bins; ibin++) {
    x[ibin - 1] = fHisto2D->GetXaxis()->GetBinCenter(ibin);
    xerr[ibin - 1] = fHisto2D->GetXaxis()->GetBinWidth(ibin) / 2;
    auto fProjHisto2D = fHisto2D->ProjectionY("", ibin, ibin);
    y[ibin - 1] = fProjHisto2D->GetMean();
    yerr[ibin - 1] = fProjHisto2D->GetMeanError();
  }
  TGraphErrors *fGraphOut = new TGraphErrors(num_of_bins, x, y, xerr, yerr);
  return fGraphOut;
}

template <typename T> TGraphErrors *MakeSigmaGraphFromTH2(T fHisto2D) {
  int num_of_bins = fHisto2D->GetXaxis()->GetLast();
  float x[num_of_bins], y[num_of_bins], xerr[num_of_bins], yerr[num_of_bins];
  for (int ibin = 1; ibin <= num_of_bins; ibin++) {
    x[ibin - 1] = fHisto2D->GetXaxis()->GetBinCenter(ibin);
    xerr[ibin - 1] = fHisto2D->GetXaxis()->GetBinWidth(ibin) / 2;
    auto fProjHisto2D = fHisto2D->ProjectionY("", ibin, ibin);
    y[ibin - 1] = fProjHisto2D->GetStdDev();
    yerr[ibin - 1] = fProjHisto2D->GetStdDevError();
  }
  TGraphErrors *fGraphOut = new TGraphErrors(num_of_bins, x, y, xerr, yerr);
  return fGraphOut;
}