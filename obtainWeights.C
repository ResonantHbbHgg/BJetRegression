{
  // weights_from_data  
  TFile* dataFile=TFile::Open("./Radion_Data_default_CSV.root");
  TTree* tree_data=(TTree*)dataFile->Get("myTrees");

  TFile* csFile=TFile::Open("../Radion_DataCS_default_CSV.root");
  TTree* tree_cs=(TTree*)csFile->Get("myTrees");
  
  // pt gamma 
  tree_data->Draw("ptPhot1:ptPhot2>>h2D_pt_data(35,20.,160.,35,20.,160.)","weight*(massggnewvtx>100 && massggnewvtx<180)*(massggnewvtx<115 || massggnewvtx>135)","colz");
  tree_cs->Draw("ptPhot1:ptPhot2>>h2D_pt_cs(35,20.,160.,35,20.,160.)","weight*(massggnewvtx>100 && massggnewvtx<180)","colz");
  float integral_2D_pt_data = h2D_pt_data->Integral();
  h2D_pt_data->Scale(1./integral_2D_pt_data);
  float integral_2D_pt_cs = h2D_pt_cs->Integral();
  h2D_pt_cs->Scale(1./integral_2D_pt_cs);
  TH2F* cs_norm_pt = h2D_pt_cs;
  h2D_pt_data->Divide(cs_norm_pt);
  h2D_pt_data->GetXaxis()->SetTitle("p_{T} Sublead Photon");
  h2D_pt_data->GetYaxis()->SetTitle("p_{T} Lead Photon");
  h2D_pt_data->SaveAs("scales_2D_pt_data_4GeVbinning.root");
}
