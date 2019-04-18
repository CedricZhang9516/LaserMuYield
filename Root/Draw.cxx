#include "/Users/zhangce/WorkArea/CZhang/CZhang.cxx"

using namespace std;

TH1D * h[5];
TH2D * hxy;
TH2D * hxxp;
TH2D * hyyp;

void Draw(){

	std::string hname = "FinalRho";
	int i = 0;
	h[i] = new TH1D(Form("%s_%d",hname.c_str(),i),Form("%s_%d",hname.c_str(),i),100,0,2);
	for(i = 1; i<5; i++){
		h[i] = new TH1D(Form("%s_%d",hname.c_str(),i),Form("%s_%d",hname.c_str(),i),2000,-1e-5,1e-5);
	}

	hxy = new TH2D("hxy","beam size XY;X(mm);Y(mm)",200,-100,100,200,-100,100);
	hxxp = new TH2D("hxxp","XXp;X(mm);X'",100,-100,100,100,-2,2);
	hyyp = new TH2D("hyyp","YYp;Y(mm);Y'",100,-100,100,100,-2,2);

	//TFile * f = new TFile("./test.root");
	//TFile * f = new TFile("./Detail_output_0.root");
	//TFile * f = new TFile("./output_0.root");

	//TFile * f = new TFile("./Detail_output_0410_0.root");
	TFile * f = new TFile("./output_0411_0.root");

	//TFile * f = new TFile("./No_Tsf/Detail_output_0.root");
	//TFile * f = new TFile("output_0410_0.root");
	//TFile * f = new TFile("./output0410_0.root");
	//TFile * f = new TFile("./test_900.root");
	TTree * t = (TTree*) f->Get("MuTree");

	//TCanvas *c2 = new TCanvas("c2","c2",800,800);
	//t->Draw( "FinalRho[][0]+FinalRho[][3]+FinalRho[][4]:FinalTime", "", "L");
	//t->Draw( "FinalRho[][4]:FinalTime", "", "");
	//t->Draw( "Rho[][0]+Rho[][3]+Rho[][4]:Time", "", "L");
	//t->Draw( "Rho[][0]:Time", "", "");


	//std::array<double,5> Rho;
	std::array<double,5> FinalRho;
	double  FinalTime;
	const int fNPoints = 201;
	double* fTime;     // (t) array
  	double (*fRho)[5];

  	double fFinalPos[3];
  	double fFinalVelocity[3];

  	fTime = new double[fNPoints]();
  	fRho = new double[fNPoints][5]();

	t->SetBranchAddress("FinalRho",&FinalRho);
	t->SetBranchAddress("FinalTime",&FinalTime);
	t->SetBranchAddress( "FinalPosition", fFinalPos );
	//t->SetBranchAddress( "FinalTime", &fFinalTime );
	t->SetBranchAddress( "FinalVelocity", fFinalVelocity );
	t->SetBranchAddress("Time",fTime);
	t->SetBranchAddress("Rho",fRho);
	//fTree->Branch( "FinalRho", fFinalRho, "FinalRho[5]/D" );


	
	TProfile * tp_3 = new TProfile("","N_tot population_2S_rho[3]; time (us); N",220,0,2.2e-6,0,1);
	TProfile * tp_4 = new TProfile("","N_tot population_ion_rho[4]; time (us); N",220,0,2.2e-6,0,1);

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
    		if(FinalRho[j]>=4.94066e-324)h[j]->Fill(FinalRho[j]);
    	}

    	for(int j = 0; j<fNPoints; j++){
    		//tp->Fill(fTime[j],fRho[j][0]);
    		//tp->Fill(fTime[j],fRho[j][0]+fRho[j][1]+fRho[j][2]+fRho[j][3]+fRho[j][4]);
    		tp_3->Fill(fTime[j],fRho[j][3]);
    		tp_4->Fill(fTime[j],fRho[j][4]);
    		//cout<<fTime[j]<<" "<<fRho[j][0]<<endl;
    	}

    	hxy->Fill(fFinalPos[0],fFinalPos[1]);
    	hxxp->Fill(fFinalPos[0],fFinalVelocity[0]/fFinalVelocity[2]);
    	hyyp->Fill(fFinalPos[1],fFinalVelocity[1]/fFinalVelocity[2]);
    }

    cout<<max1<<" "<<min1<<endl;
    cout<<max2<<" "<<min2<<endl;
	
	TCanvas *c = new TCanvas("c","c",2000,400);
  	c->Divide(5,1);
  	//c->Divide(3,1);
  	for (int i = 1; i<6; i++){c->cd(i);SetstyleHist1(h[i-1]);h[i-1]->Draw();}
  	//c->cd(1);h[0]->Draw();
  	//c->cd(2);h[3]->Draw();
  	//c->cd(3);h[4]->Draw();
  	TCanvas *c3 = new TCanvas("c3","c3",800,400);
  	c3->Divide(2,1);
  	c3->cd(1);
  	tp_3->Draw();
  	c3->cd(2);
  	tp_4->Draw();

  	TCanvas *c4 = new TCanvas("c4","c4",1200,400);
  	c4->Divide(3,1);
  	c4->cd(1);
  	hxy->Draw("colz");

  	c4->cd(2);
  	hxxp->Draw("colz");

  	c4->cd(3);
  	hyyp->Draw("colz");



}



