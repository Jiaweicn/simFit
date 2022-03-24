#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "Fit/LogLikelihoodFCN.h"
#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include "Math/FitMethodFunction.h"
#include "Fit/BasicFCN.h"
#include "Fit/FcnAdapter.h"
#include "Fit/FitConfig.h"
#include "Fit/FitResult.h"
#include "Math/Error.h"
#include <memory>
#include "Math/IParamFunction.h"
#include "Math/MultiDimParamFunctionAdapter.h"
#include <fstream>
#include <iostream>

clock_t tstart, tend;//output run time of code

Int_t i, j;

//array of all need-to-fit histograms storing in root file
string hlist[] = {"h1011","h1012","h1021","h1022","h1023","h1024","h1031","h1032","h1033","h1034","h1041","h1042","h1043","h1044","h1051","h1052","h1053","h1054","h1061","h1062","h1063","h1064","h1071","h1072","h1073","h1074"};
const Int_t nhists=sizeof(hlist)/sizeof(hlist[0]);

//degree values in C.M. frame of corresponding histograms
Double_t theta[]={0.742,1.025,2.368,2.815,3.275,3.499,3.734,3.97,4.429,4.888,5.241,5.712,6.172,6.642,6.996,7.455,7.925,8.396,8.749,9.22,9.69,10.16,10.502,10.972,11.442,11.912};

char rootinput[20]="dsing_h_shift.root";
TFile *fin= TFile::Open(rootinput);//input histos which need to be fitted
char rootoutname[30]="simFitResult.root";
TFile *rootout = TFile::Open(rootoutname,"UPDATE");//output fitted histos into a root file

//Double_t xmin=13.2, xmax=16.27;
//Double_t xmin=16.3, xmax=19.;
//Double_t xmin=19., xmax=20.63;
Double_t xmin=18., xmax=18.3;

//array of peaks' location
Double_t pks[]={
18.15
	/*6.44,
	6.6,
	6.71,
	6.99,
	7.15,
	7.36,
	7.63,
	8.36,
	8.44,
	9.02,
	9.15,
	9.31,
	9.49,
	9.83,
	10.03*/

	/*10.35,
	10.6,
	10.68,
	10.91,
	11.01,
	11.18,
	11.39,
	11.52,
	11.73,
	11.86,
	11.99,
	12.26*/
	
	//12.5,
	//12.79,
	//13.11,
	/*13.34,
	13.4,
	13.72,
	13.86,
	14.10,
	14.25,
	14.47,
	14.72,
	14.81,
	15.03,
	15.15,
	15.31,
	15.46,
	15.75,
	16.11*/
	
	/*16.49,
	16.67,
	16.82,
	16.92,
	17.33,
	17.56,
	17.7,
	18.15,
	18.56,
	18.66*/

	/*19.4,
	19.5,
	19.66,
	19.75,
	19.85,
	20.11,
	20.22,
	20.4*/

};
const Int_t npks=sizeof(pks)/sizeof(Double_t);


const Int_t nbak = 3;//number of parameters of background function
const Int_t nind = nbak+npks;//number of NOT-shared parameters in each histogram
const Int_t nshr = 2*npks;//number of shared parameters in each histogram
const Int_t n = nind+nshr;//total number of parameters in each histogram


//setting shared parameters
Int_t par_h[nhists][n];//={{0,1,2,3,4,15,16,17,18}, {5,6,7,8,9,15,16,17,18} ,{10,11,12,13,14,15,16,17,18}};

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Lorenzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
	return (0.5*par[0]*par[2]/TMath::Pi()) /TMath::Max(1.e-10,(x[0]-par[1])*(x[0]-par[1])+0.25*par[2]*par[2]);
}

// normalized Gaus peak
Double_t gausPeak(Double_t *x, Double_t *par) {
	return TMath::Gaus(x[0],par[1],par[2], 1)*par[0];
}

const Double_t sigma=0.02;//in MeV
Double_t voigtPeak(Double_t *x, Double_t *par){
	Double_t lg=par[2];// lg is FWHM of Lorentzian Distribution
	return TMath::Voigt(x[0]-par[1],sigma,lg,2)*par[0];//Voigt itself is normalized
}// sigma is FWHM of Gaussian Distribution, fixed for spread of beam.

Double_t peaksFunction(Double_t *x, Double_t *par) {
	Double_t fitval = 0;
	for(i=0;i<npks;i++) fitval=fitval+gausPeak(x,&par[3*i+0]);
	return fitval;
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
	Double_t fitval = background(x,par);
	fitval=fitval+peaksFunction(x,&par[3]);
	return fitval;
}


// Create the GlobalChi2 structure
struct Global {

	//ROOT::Math::IMultiGenFunction *fChi2[nhists];
	std::vector<const ROOT::Math::IMultiGenFunction *> fChi2;

	Global(std::vector<ROOT::Fit::Chi2Function *> &func) { 
		for(int k=0, size=func.size(); k<size; k++) fChi2.push_back(func[k]);
	}
	
	// setting parameters for each histogram //
	Double_t operator() (const Double_t *par) const {
		Double_t p[nhists][n];
		for(int k=0, size=fChi2.size(); k<size; k++) {
			for(int m=0; m<n; m++) p[k][m] = par[par_h[k][m]];
		}
		Double_t ret = 0;
		for(int k=0, size=fChi2.size(); k<size; k++) {
			ret += (*fChi2[k])(p[k]);
		}
		return ret;
   }
};

void currenttime(){
   time_t now = time(0);
   char* dt = ctime(&now);
   cout << "Current time: " << dt << endl;
}


void plotFit(TH1F *hists[]){//plot fitting results

	TCanvas * chist = new TCanvas("cSimfit","Simultaneous fit of multi-histograms",10,10,2500,1500);
	if(nhists%4 ==0) chist->Divide(4,nhists/4);
	if(nhists%4 !=0) chist->Divide(4,1+nhists/4);
	gStyle->SetOptFit(0000);
	gStyle->SetOptStat(0);
	for(i=0; i<nhists; i++){
		chist->cd(i+1);
		hists[i]->Draw();
	}

	chist->Write("",TObject::kOverwrite);
}

//print out all parameters to screen
void par2Screen(TF1 *fitfs[]){
	cout<<"-------++++++++++Chi2 fitting result+++++++++--------"<<endl;
	cout<<"Locations and Widthes of all Peaks"<<endl;
	//print out all peaks' parameters
	for(i=0; i<npks; i++){
		j=nbak+3*i+1;
		cout<<"Location: "<<fitfs[0]->GetParameter(j)<<" +/- "<<fitfs[0]->GetParError(j)<<"  Width: "<<fitfs[0]->GetParameter(j+1)<<" +/- "<<fitfs[0]->GetParError(j+1)<<endl;
	}

	cout<<"Strengthes and other individual parameters"<<endl;
	//print out strengthes and other parameters
	for(i=0; i<nhists; i++){
		cout<<"bkg: ";
		for(j=0; j<nbak; j++){
		cout<<fitfs[i]->GetParameter(j)<<" +/- "<<fitfs[i]->GetParError(j)<<"  ";
		}
		cout<<"strength: ";
		for(j=0; j<npks; j++){
		cout<<fitfs[i]->GetParameter(nbak+3*j)<<" +/- "<<fitfs[i]->GetParError(nbak+3*j)<<"  ";
		}
		cout<<""<<endl;
	}
}
//print out all parameters to a file
void par2File(char txtname[30],TF1 *fitfs[]){
	ofstream fout;
	fout.open(txtname);

	cout<<"Output paras txt file: "<<txtname<<endl;
	fout<<"Fit_range: "<<xmin<<" "<<xmax<<endl;
	fout<<"Locations and Widthes of all Peaks"<<endl;
	//print out all peaks' parameters
	for(i=0; i<npks; i++){
		j=nbak+3*i+1;
		fout<<"Location: "<<fitfs[0]->GetParameter(j)<<" +/- "<<fitfs[0]->GetParError(j)<<"  Width: "<<fitfs[0]->GetParameter(j+1)<<" +/- "<<fitfs[0]->GetParError(j+1)<<endl;
	}

	fout<<"Strengthes and other individual parameters"<<endl;
	//print out strengthes and other parameters
	for(i=0; i<nhists; i++){
		fout<<"bkg: ";
		for(j=0; j<nbak; j++){
		fout<<fitfs[i]->GetParameter(j)<<" +/- "<<fitfs[i]->GetParError(j)<<"  ";
		}
		fout<<"strength: ";
		for(j=0; j<npks; j++){
		fout<<fitfs[i]->GetParameter(nbak+3*j)<<" +/- "<<fitfs[i]->GetParError(nbak+3*j)<<"  ";
		}
		fout<<""<<endl;
	}

}


void angDis(TF1 *fitfuncs[]){//plot angular distribution from fit result
	TCanvas *cgre = new TCanvas("angDist0","Ang-Dist from Chi2 fitting",10,10,1500,1000);
	TGraphErrors *gre[npks];
	if(npks%4 ==0) cgre->Divide(4,npks/4);
	if(npks%4 !=0) cgre->Divide(4,1+npks/4);
	Double_t num;
	char grename[40];
	char gretitle[30];
	for(i=0; i<npks; i++){
		gre[i]=new TGraphErrors(nhists);

		for(j=0; j<nhists; j++){
			gre[i]->SetPoint(j,theta[j],fitfuncs[j]->GetParameter(nbak+3*i));
			gre[i]->SetPointError(j,0.,fitfuncs[j]->GetParError(nbak+3*i));
		}
		cgre->cd(i+1)->SetLogy();
		
		num=fitfuncs[0]->GetParameter(nbak+3*i+1);
		sprintf(grename,"angDist0_%.0fdot%.0f", floor(num), 1000*(num-floor(num)));
		sprintf(gretitle,"%.3f MeV peak",num);
		gre[i]->SetName(grename);
		gre[i]->SetTitle(gretitle);
		
		gre[i]->SetMarkerColor(4);
		gre[i]->SetMarkerStyle(21);
		gre[i]->Draw("ALP");
		//gre[i]->Write("",TObject::kOverwrite);
	}
	cgre->Write("",TObject::kOverwrite);
}

void simFit(){

	tstart=clock();

	currenttime();
	
	if (fin == nullptr) {
		cout<<"input root file wrong! abort!"<<endl;
		return 1;
	}
	
	//ROOT::EnableImplicitMT(8);//run code in multiple cores
	cout<<"Input root file: "<<rootinput<<endl;
	cout<<"Fitting range: "<<xmin<<" to "<<xmax<<endl;
	cout<<"Number of peaks: "<<npks<<endl;
	cout<<"Number of histograms: "<<nhists<<endl;
	cout<<"Output root file: "<<rootoutname<<endl;
	
	for(i=0; i<n; i++) par_h[0][i]=i;
	for(i=1; i<nhists; i++) {
		for(j=0; j<nbak; j++) par_h[i][j]=n+nind*(i-1)+j;
		for(j=0; j<npks; j++){
			par_h[i][nbak+3*j]=n+nind*(i-1)+nbak+j;
			par_h[i][nbak+3*j+1]=par_h[0][nbak+3*j+1];
			par_h[i][nbak+3*j+2]=par_h[0][nbak+3*j+2];
		}
	}

	TH1F *hists[nhists];//an array of histograms
	TF1 *fitfcns[nhists];
	ROOT::Math::WrappedMultiTF1 *wfs[nhists];
	ROOT::Fit::DataOptions opt;
	ROOT::Fit::DataRange range[nhists];
	ROOT::Fit::BinData *data[nhists];
	//ROOT::Fit::UnBinData *data[nhists];
	ROOT::Fit::Chi2Function *chi2[nhists];
	std::vector<ROOT::Fit::Chi2Function *> chi2_vec;
	//ROOT::Fit::PoissonLLFunction ll[nhists];//Likelyhood fit method

	char hname[10];//char variable used to load histogram in every step
	char fname[10];//fit function name
	for(i=0 ; i<nhists; i++){

		strcpy(hname, hlist[i].c_str());
		fin->GetObject(hname, hists[i]);//get the histogram from root file
		hists[i]->GetXaxis()->SetRangeUser(xmin,xmax);//set needed range of X axis of histogram
		sprintf(fname,"f_%s",hname);//function name of fitting in every histogram
		fitfcns[i] = new TF1(fname, fitFunction, xmin, xmax, n);
		wfs[i] = new ROOT::Math::WrappedMultiTF1(*fitfcns[i],fitfcns[i]->GetNdim());

		range[i].SetRange(xmin, xmax);
		data[i] = new ROOT::Fit::BinData(opt,range[i]);
		//data[i] = new ROOT::Fit::UnBinData(opt,range[i]);
		ROOT::Fit::FillData(*data[i], hists[i]);

		chi2[i] = new ROOT::Fit::Chi2Function(*data[i], *wfs[i]);
		chi2_vec.push_back(chi2[i]);
		//ll[i] = new ROOT::Fit::PoissonLLFunction(*data[i], *wfs[i]);
	}
	// perform now global fit
	Global global(chi2_vec);

	ROOT::Fit::Fitter fitter;

	//setting parameters
	Int_t Npar = nhists*3 + 2*npks + nhists*npks;
	Double_t par0[Npar];

	// create before the parameter settings in order to fix or set range on them
	fitter.Config().SetParamsSettings(Npar,par0);

	//set the parameters of all histograms (background & strength of peaks)
	char pname[20];//name of parameters

	//setting not-shared parameters for 1st histogram
	sprintf(pname, "const_%d",1);//name of constant term of quadratic background
	fitter.Config().ParSettings(0).Set(pname, 10., 0.1, 0., 1.e3);
	sprintf(pname, "linear_%d",1);
	fitter.Config().ParSettings(1).Set(pname,0);
	fitter.Config().ParSettings(1).Fix();
	sprintf(pname, "quad_%d",1);
	fitter.Config().ParSettings(2).Set(pname,0);
	fitter.Config().ParSettings(2).Fix();
	for(i=0; i<npks; i++){//strength parameters for 1st histogram
		sprintf(pname, "strength1-%d",i+1);
		fitter.Config().ParSettings(nbak+3*i).Set(pname, 30., 0.1, 0., 1.e3);
	}
	// setting shared parameters of all peaks 
	for(i=0; i<npks; i++){
		sprintf(pname, "peak_%d",i+1);
		fitter.Config().ParSettings(nbak+3*i+1).Set(pname, pks[i], 0.001, pks[i]-0.05 ,pks[i]+0.05);
		sprintf(pname, "width_%d",i+1);
		fitter.Config().ParSettings(nbak+3*i+2).Set(pname, 0.05, 0.001, 1e-3,5.e-1);
	}
	
	//setting not-shared parameters for other histograms
	for(i=1; i<nhists; i++){
		sprintf(pname, "const_%d",i+1);//name of constant term of quadratic background
		fitter.Config().ParSettings(n+nind*(i-1)+0).Set(pname, 10., 0.1, 0., 1.e3);//set the constant term of background
		sprintf(pname, "linear_%d",i+1);//name of linear term of quadratic background
		fitter.Config().ParSettings(n+nind*(i-1)+1).Set(pname,0);
		fitter.Config().ParSettings(n+nind*(i-1)+1).Fix();//fix the linear term ofbackground to 0
		sprintf(pname, "quad_%d",i+1);//name of quadratic term of quadratic background
		fitter.Config().ParSettings(n+nind*(i-1)+2).Set(pname,0);
		fitter.Config().ParSettings(n+nind*(i-1)+2).Fix();//fix the quadratic term of background to 0
		for(j=0; j<npks; j++){
			sprintf(pname, "strength%d-%d",i+1, j+1);
			fitter.Config().ParSettings(n+nind*(i-1)+nbak+j).Set(pname, 30., 0.1, 0., 1.e3);
		}
	}
 
	fitter.Config().MinimizerOptions().SetPrintLevel(0);
	fitter.Config().SetMinimizer("Minuit","Migrad");

	// fit FCN function directly
	// (specify optionally data size and flag to indicate that is a chi2 fit)
	Int_t size=0;
	for(i=0; i<nhists; i++){
		size=size+data[i]->Size();
	}
	fitter.FitFCN(Npar,global, 0, size, true);
	ROOT::Fit::FitResult result = fitter.Result();
	
	//result.Print(std::cout);

	for(i=0; i<nhists; i++){
		fitfcns[i]->SetFitResult(result, par_h[i]);
		fitfcns[i]->SetLineColor(kRed);
		fitfcns[i]->SetNpx(3000);
		hists[i]->GetListOfFunctions()->Add(fitfcns[i]);
		//hists[i]->Write("",TObject::kOverwrite);
	}

	//plot the fitting results
	plotFit(hists);

	//print out parameters
	char txtname[30];
	//sprintf(txtname,"simFitResult%.2fto%.2f.txt",xmin,xmax);
	sprintf(txtname,"simFitResult.txt");
	par2File(txtname, fitfcns);
	//par2Screen(fitfcns);

	//draw the angular distribution
	//angDis(fitfcns);
	

	char listname[20]="histolist";
	TList * l= new TList();
	for(i=0; i<nhists; i++)
	{
		l->Add(hists[i]);
	}
	//l->Write(listname,TObject::kOverwrite);
	l->Write(listname, TObject::kSingleKey);
	rootout->Close();

	Double_t chisqr=fitfcns[0]->GetChisquare();
	Double_t ndf=fitfcns[0]->GetNDF();
	cout<<"-----+++ Goodness of the fit +++-----"<<endl;
	cout<<"Probability: "<<fitfcns[0]->GetProb()<<endl;
	cout<<"Chi2/NDF: "<<chisqr<<"/"<<ndf<<"="<<chisqr/ndf<<endl;
	cout<<"------++++*****$$END$$*****++++------"<<endl;
	currenttime();
	cout<<endl;

	tend=clock();//in seconds

	Double_t time_used =((double)(tend-tstart))/CLOCKS_PER_SEC;
	cout<<"run time of this code: "<<time_used/60<<" minutes | "<<time_used/60/60<<" hours"<<"\n"<<endl;
}
