//used to draw raw and corrected angular distribution from files
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include <sstream>
//angles are in C.M. frame

Double_t thetaCM[]={//scattering angle of all histograms in C.M. frame
	0.742,
	1.025,
	2.368,
	2.815,
	3.275,
	3.499,
	3.734,
	3.97,
	4.429,
	4.888,
	5.241,
	5.712,
	6.172,
	6.642,
	6.996,
	7.455,
	7.925,
	8.396,
	8.749,
	9.22,
	9.69,
	10.16,
	10.502,
	10.972,
	11.442,
	11.912
};
const Int_t nthetaCM=sizeof(thetaCM)/sizeof(thetaCM[0]);//number of thetaCM equals to number of histograms


Int_t i,j;
Double_t num;
vector<vector<Double_t>> err;//load errors of all strengthes of all peaks
vector<Double_t> pks;
vector<Double_t> strth0;//load strengthes of all peaks
vector<Double_t> err0;//load errors of all strengthes of all peaks
vector<vector<Double_t>> strth;//load strengthes of all peaks
Double_t xmin, xmax;
char file[20]="simFitResult";
char rootname[30];
char xtitle[20]="#theta_{c.m.}";//set the title of X axis
char ytitle[20]="XSection";
const Int_t nbaks=3;//number of parameters in background part of fitting function
Int_t npks;
Int_t nhists;

string theorylist[]={"j0.txt","j1.txt","j2.txt","j3.txt","j4.txt",};
const Int_t ntheory=sizeof(theorylist)/sizeof(theorylist[0]);
TGraph *theoryGr[ntheory];//array of theoretical angular distrbution
void plottheory(){
	const Double_t scale=0.05;
	TMultiGraph *mg = new TMultiGraph();
	char dir[20]="../theoryAngDist/";//directory of theoretical angular distrbution
	ifstream theoryfile;
	string line;//used to read a whole line out of the file
	char filename0[10];
	char filename[50];
	Double_t x, y;
	for(i=0; i<ntheory; i++){
		strcpy(filename0,theorylist[i].c_str());
		sprintf(filename, "%s%s", dir, filename0);
		theoryfile.open(filename);
		theoryGr[i]=new TGraph();
		j=0;
		while(getline(theoryfile, line)){//get a line of data from file to variable fin
			stringstream ss(line);
			ss >> x;
			ss >> y;
			theoryGr[i]->SetPoint(j, x, scale*y);
			j++;
		}
		theoryGr[i]->SetMarkerColor(i+1);
		theoryGr[i]->SetLineColor(i+1);
		theoryGr[i]->SetMarkerStyle(1);
		//theoryGr[i]->SetMarkerSize(0.3);
		theoryfile.close();
	}
	cout<<"Theoretical angular distrbutions plotted!!!"<<endl;
}


void aboutfile(){//determine the number of peaks in the fitting file
	
	ifstream fin(file);
	string line;//used to read a whole line out of the file
	string unwanted;//used to load string in this line
	i=j=0;
	nhists=0;
	while(getline(fin, line)){//get a line of data from file to variable fin
		stringstream ss(line);
		ss >> unwanted;//used to take string "bkg:" and "Location:" and "strength"
		if(unwanted=="Fit_range:"){
			ss>>xmin;
			ss>>xmax;
		}
		if(unwanted=="Location:"){
			ss>>num;
			pks.push_back(num);
		}
		npks=pks.size();
		if(unwanted=="bkg:"){
			for(j=0;j<3*nbaks;j++) ss>>unwanted;
			ss>>unwanted;
			if(unwanted=="strength:"){
				for(j=0;j<npks;j++){
					ss>>num;
					strth0.push_back(num);
					ss>>unwanted;
					ss>>num;
					err0.push_back(num);
				}
				strth.push_back(strth0);
				err.push_back(err0);
			}
		}
		strth0.clear();
		err0.clear();
		i++;
	}
	nhists=strth.size();
	if(nthetaCM!=nhists) cout<<"number of histos != number of angles"<<endl;
	fin.close();
	cout<<"number of peaks: "<<npks<<endl;
	cout<<"number of histograms: "<<nhists<<endl;
}


void rawAngDist(){//angular distribution of raw fit strengthes of all peaks
	TCanvas *cgre = new TCanvas("angularDistribution","Raw Fit Angular Distribution",10,10,1500,1000);
	TMultiGraph *mg[npks];
	TGraphErrors *gre[npks];
	if(npks%4 ==0) cgre->Divide(4,npks/4);
	if(npks%4 !=0) cgre->Divide(4,1+npks/4);
	char title[30];
	for(i=0; i<npks; i++){
		mg[i]=new TMultiGraph();
		gre[i]=new TGraphErrors(nhists);
		for(j=0; j<nhists; j++){
			gre[i]->SetPoint(j,thetaCM[j],strth[j][i]);
			gre[i]->SetPointError(j,0.,err[j][i]);
		}
		cgre->cd(i+1)->SetLogy();
		sprintf(title,"%.3f MeV peak",pks[i]);
		gre[i]->SetLineColor(6);
		gre[i]->SetMarkerColor(6);
		gre[i]->SetMarkerStyle(29);
		//gre[i]->Draw("ALP");
		mg[i]->Add(gre[i]);
		for(j=0; j<ntheory; j++){
			mg[i]->Add(theoryGr[j]);
		}
		mg[i]->SetTitle(title);
		mg[i]->Draw("ALP");
		mg[i]->GetXaxis()->SetTitle(xtitle);
		mg[i]->GetYaxis()->SetTitle(ytitle);
	}
}

void angDist6(){//give the angular distribution corrected by Kawabata 6.433MeV angular distribution data
	Double_t thetaCM_kawa[]={
	0.742,
	//1.025,//
	2.368,
	2.815,
	3.275,
	3.499,
	//3.734,//
	3.97,
	4.429,
	4.888,
	5.241,
	5.712,
	6.172,
	6.642,
	6.996,
	7.455,
	7.925,
	8.396,
	8.749,
	9.22,
	9.69,
	10.16,
	10.502,
	10.972,
	11.442,
	11.912
	};
	Double_t kawa6dot433[]={//the 2ed and 7th data point of 6.43 0+ state are abandoned in Kawabata paper
	67,
	//0.1,//
	9,
	2.8,
	0.9,
	1.6,
	//0.1,//
	3.5,
	6.5,
	7.0,
	7.1,
	5.2,
	5.1,
	2.,
	0.9,
	0.36,
	0.57,
	1.1,
	1.5,
	1.6,
	1.5,
	0.95,
	0.7,
	0.41,
	0.42,
	0.22
	};
	Double_t my6dot433[]={
	396.27,
	46.9855,
	13.5793,
	4.76335,
	5.66937,
	11.4169,
	21.0286,
	25.3607,
	47.37,
	41.0301,
	30.5758,
	18.1782,
	15.0834,
	5.26564,
	5.5477,
	11.0855,
	59.7498,
	61.2413,
	54.9194,
	38.3618,
	24.494,
	14.7049,
	23.2487,
	11.6861
	};
	const Int_t nkawa=sizeof(kawa6dot433)/sizeof(kawa6dot433[0]);
	TCanvas *ckawa = new TCanvas("angDist6","Corrected Angular Distribution(6.433MeV)",10,10,1500,1000);
	TMultiGraph *mg[npks];
	TGraphErrors *gre[npks];
	//TGraph *gr[npks];
	if(npks%4 ==0) ckawa->Divide(4,npks/4);
	if(npks%4 !=0) ckawa->Divide(4,1+npks/4);
	char title[30];
	for(i=0; i<npks; i++){
		mg[i]=new TMultiGraph();
		gre[i]=new TGraphErrors(nkawa);
		for(j=0; j<nkawa; j++){
			num=kawa6dot433[j]/my6dot433[j];//nbak+6 is the peak position of 6.43
			if(j>=6){
				gre[i]->SetPoint(j,thetaCM_kawa[j],num*strth[j+2][i]);
				gre[i]->SetPointError(j,0.,num*err[j+2][i]);
			}
			else if(j>=1 && j<6){
				gre[i]->SetPoint(j,thetaCM_kawa[j],num*strth[j+1][i]);
				gre[i]->SetPointError(j,0.,num*err[j+1][i]);
			}
			else{
				gre[i]->SetPoint(j,thetaCM_kawa[j],num*strth[j][i]);
				gre[i]->SetPointError(j,0.,num*err[j][i]);
			}
		}
		ckawa->cd(i+1)->SetLogy();
		sprintf(title,"%.3f MeV peak",pks[i]);
		gre[i]->SetLineColor(6);
		gre[i]->SetMarkerColor(6);
		gre[i]->SetMarkerStyle(29);
		//gre[i]->Draw("ALP");
		mg[i]->Add(gre[i]);
		for(j=0; j<ntheory; j++){
			mg[i]->Add(theoryGr[j]);
		}
		mg[i]->SetTitle(title);
		mg[i]->Draw("ALP");
		mg[i]->GetXaxis()->SetTitle(xtitle);
		mg[i]->GetYaxis()->SetTitle(ytitle);
	}
	TFile *rootout = TFile::Open(rootname,"UPDATE");//output into a root file
	ckawa->Write("",TObject::kOverwrite);
	rootout->Close();
}


void angDist9(){//give the angular distribution normalized by Kawabata 9.305MeV angular distribution data
	Double_t kawa9dot305[]={// 9.305MeV peak angular distribution data from Kawabata paper
	18.59,
	16.89,
	3.07,
	0.91,
	0.47,
	0.78,
	1.29,
	1.90,
	2.77,
	2.75,
	2.20,
	1.94,
	0.90,
	0.53,
	0.52,
	0.62,
	0.87,
	1.13,
	1.45,
	1.31,
	1.11,
	0.76,
	0.59,
	0.44,
	0.35,
	0.41 
	};
	Double_t my9dot305[]={
	130.98,
	105.56,
	15.61,
	3.75,
	1.77,
	5.40,
	2.31,
	5.50,
	8.17,
	7.86,
	11.11,
	9.54,
	4.51,
	2.37,
	3.96,
	5.12,
	7.47,
	11.20,
	51.97,
	47.35,
	38.83,
	24.91,
	17.87,
	12.97,
	9.44,
	11.65,
	};
	const Int_t nkawa=sizeof(kawa9dot305)/sizeof(kawa9dot305[0]);
	TCanvas *ckawa = new TCanvas("angDist9","Corrected Angular Distribution(9.305MeV)",10,10,1500,1000);
	TMultiGraph *mg[npks];
	TGraphErrors *gre[npks];
	if(npks%4 ==0) ckawa->Divide(4,npks/4);
	if(npks%4 !=0) ckawa->Divide(4,1+npks/4);
	char title[30];
	for(i=0; i<npks; i++){
		mg[i]=new TMultiGraph();
		gre[i]=new TGraphErrors(nkawa);
		for(j=0; j<nkawa; j++){
			num=kawa9dot305[j]/my9dot305[j];
			gre[i]->SetPoint(j,thetaCM[j],num*strth[j][i]);
			gre[i]->SetPointError(j,0.,num*err[j][i]);
		}
		ckawa->cd(i+1)->SetLogy();
		sprintf(title,"%.3f MeV peak",pks[i]);
		gre[i]->SetLineColor(6);
		gre[i]->SetMarkerColor(6);
		gre[i]->SetMarkerStyle(29);
		//gre[i]->Draw("ALP");
		mg[i]->Add(gre[i]);
		for(j=0; j<ntheory; j++){
			mg[i]->Add(theoryGr[j]);
		}
		mg[i]->SetTitle(title);
		mg[i]->Draw("ALP");
		mg[i]->GetXaxis()->SetTitle(xtitle);
		mg[i]->GetYaxis()->SetTitle(ytitle);
	}
	TFile *rootout = TFile::Open(rootname,"UPDATE");//output into a root file
	ckawa->Write("",TObject::kOverwrite);
	rootout->Close();
}

void plotAngDist(){
	//sprintf(file,"%s%.2fto%.2f.txt",file,xmin,xmax);//specify the input file name
	sprintf(rootname,"%s.root",file);//specify the input file name
	sprintf(file,"%s.txt",file);//specify the input file name
	cout<<"input txt file: "<<file<<endl;
	cout<<"output root file: "<<rootname<<endl;

	aboutfile();//to get the info of numbers of peaks and histograms
	cout<<"fit range of file: "<<xmin<<" to "<<xmax<<endl;
	
	plottheory();//fill theoretical angular distrbutions
	
	//rawAngDist();
	angDist6();
	angDist9();
}
