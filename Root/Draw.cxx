#include "/Users/zhangce/WorkArea/CZhang/CZhang.cxx"

using namespace std;

TH1D * h[5];

void Draw(){

	std::string hname = "Rho";
	int i = 0;
	h[i] = new TH1D(Form("%s_%d",hname.c_str(),i),Form("%s_%d",hname.c_str(),i),100,0,2);
	for(i = 1; i<5; i++){
		h[i] = new TH1D(Form("%s_%d",hname.c_str(),i),Form("%s_%d",hname.c_str(),i),200,-1e-5,1e-5);
	}

	//TFile * f = new TFile("./test.root");
	//TFile * f = new TFile("./Detail_output_0.root");
	//TFile * f = new TFile("./output_0.root");
	TFile * f = new TFile("./Detail_output_0410_0.root");
	//TFile * f = new TFile("./No_Tsf/Detail_output_0.root");
	//TFile * f = new TFile("output_0410_0.root");
	//TFile * f = new TFile("./output0410_0.root");
	//TFile * f = new TFile("./test_900.root");
	TTree * t = (TTree*) f->Get("MuTree");

	TCanvas *c2 = new TCanvas("c2","c2",800,800);
	//t->Draw( "FinalRho[][0]+FinalRho[][3]+FinalRho[][4]:FinalTime", "", "L");
	//t->Draw( "FinalRho[][4]:FinalTime", "", "");
	//t->Draw( "Rho[][0]+Rho[][3]+Rho[][4]:Time", "", "L");
	t->Draw( "Rho[][0]:Time", "", "");


	//std::array<double,5> Rho;
	std::array<double,5> FinalRho;
	double  FinalTime;
	const int fNPoints = 201;
	double* fTime;     // (t) array
  	double (*fRho)[5];

  	fTime = new double[fNPoints]();
  	fRho = new double[fNPoints][5]();

	t->SetBranchAddress("FinalRho",&FinalRho);
	t->SetBranchAddress("FinalTime",&FinalTime);
	t->SetBranchAddress("Time",fTime);
	t->SetBranchAddress("Rho",fRho);
	//fTree->Branch( "FinalRho", fFinalRho, "FinalRho[5]/D" );


	
	TProfile * tp = new TProfile("","N_tot population; time (us); N",100,0,2.2e-6,0,1);

	double min1 = 1, max1 = 0;
	double min2 = 0, max2 = 0;
	const int nEntries = t->GetEntries();
  	for( int entry=0; entry<nEntries; entry++ ){
  		t->GetEntry( entry );
  		//cout<<Rho[3]<<Rho[4]<<" "<<FinalTime<<endl;
  		if(FinalRho[4]>max1)max1=FinalRho[4];
  		if(FinalRho[4]>0 && FinalRho[4]<min1)min1=FinalRho[4];
  		if(FinalRho[3]>max2)max2=FinalRho[4];
  		if(FinalRho[3]>0 && FinalRho[3]<min2)min2=FinalRho[4];
  		//cout<<Rho[0]<<" "<<Rho[1]<<" "<<Rho[2]<<" "<<Rho[3]<<" "<<Rho[4]<<endl;
  		for(int j = 0; j<5; j++){
    		h[j]->Fill(FinalRho[j]);
    	}

    	for(int j = 0; j<fNPoints; j++){
    		//tp->Fill(fTime[j],fRho[j][0]);
    		//tp->Fill(fTime[j],fRho[j][0]+fRho[j][1]+fRho[j][2]+fRho[j][3]+fRho[j][4]);
    		tp->Fill(fTime[j],fRho[j][3]);
    		//cout<<fTime[j]<<" "<<fRho[j][0]<<endl;
    	}
    }

    cout<<max1<<" "<<min1<<endl;
    cout<<max2<<" "<<min2<<endl;
	
	TCanvas *c = new TCanvas("c","c",2000,400);
  	c->Divide(5,1);
  	//c->Divide(3,1);
  	for (int i = 1; i<6; i++){c->cd(i);h[i-1]->Draw();}
  	//c->cd(1);h[0]->Draw();
  	//c->cd(2);h[3]->Draw();
  	//c->cd(3);h[4]->Draw();
  	TCanvas *c3 = new TCanvas("c3","c3",2000,400);
  	tp->Draw();




}



