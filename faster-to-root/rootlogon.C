

{


  
  	TFile *_file0 = TFile::Open("run.root");

	TCanvas c1("c1", "c1");
  	DataTree->Draw("Si>>h(5000,50000,60000)");
 
  	TCanvas c2("c2", "c2");
 	DataTree->Draw("TAC");
  	
  	

  	TCanvas c3("c3", "c3");
 	DataTree->Draw("Si:TAC>>h2(3000,0,700000,3000,0,70000)","","col");



}
