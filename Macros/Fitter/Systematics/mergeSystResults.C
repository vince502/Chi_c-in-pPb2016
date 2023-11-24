#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "../../ChiFitterInit.h"

//Merge the fits from different files (with the same naming scheme filename_#.root) to one file

void CopyDir(TDirectory *source) {
	//copy all objects and subdirs of directory source as a subdir of the current directory
	source->ls();
	TDirectory *savdir = gDirectory;
	TDirectory *adir = savdir->mkdir(source->GetName());
	adir->cd();
	//loop on all entries of this directory
	TKey *key;
	TIter nextkey(source->GetListOfKeys());
	while ((key = (TKey*)nextkey())) {
		const char *classname = key->GetClassName();
		TClass *cl = gROOT->GetClass(classname);
		if (!cl) continue;
		if (cl->InheritsFrom(TDirectory::Class())) {
			source->cd(key->GetName());
			TDirectory *subdir = gDirectory;
			adir->cd();
			CopyDir(subdir);
			adir->cd();
		}
		else if (cl->InheritsFrom(TTree::Class())) {
			TTree *T = (TTree*)source->Get(key->GetName());
			adir->cd();
			TTree *newT = T->CloneTree(-1, "fast");
			newT->Write();
		}
		else {
			source->cd();
			TObject *obj = key->ReadObj();
			adir->cd();
			obj->Write();
			delete obj;
		}
	}
	adir->SaveSelf(kTRUE);
	savdir->cd();
}

void CopyDirFile(TDirectory *source) {
	//This is here in order to NOT replicate the file name on the top

	source->ls();
	TDirectory *adir = gDirectory;
	adir->cd();
	//loop on all entries of this directory
	TKey *key;
	TIter nextkey(source->GetListOfKeys());
	while ((key = (TKey*)nextkey())) {
		const char *classname = key->GetClassName();
		TClass *cl = gROOT->GetClass(classname);
		if (!cl) continue;
		if (cl->InheritsFrom(TDirectory::Class())) {
			source->cd(key->GetName());
			TDirectory *subdir = gDirectory;
			adir->cd();
			CopyDir(subdir);
			adir->cd();
		}
		else if (cl->InheritsFrom(TTree::Class())) {
			TTree *T = (TTree*)source->Get(key->GetName());
			adir->cd();
			TTree *newT = T->CloneTree(-1, "fast");
			newT->Write();
		}
		else {
			source->cd();
			TObject *obj = key->ReadObj();
			adir->cd();
			obj->Write();
			delete obj;
		}
	}
	adir->SaveSelf(kTRUE);
}

void CopyFile(const char *fname) {
	//Copy all objects and subdirs of file fname to the subdirectory of the current directory
	TDirectory *target = gDirectory;
	TFile *f = TFile::Open(fname);
	if (!f || f->IsZombie()) {
		printf("Cannot copy file: %s\n", fname);
		target->cd();
		return;
	}
	//f->cd(TreeNameA);
	TDirectory *dirToCopy = gDirectory;
	target->cd();
	CopyDirFile(dirToCopy);
	f->Close();
	target->cd();
}


//void mergeSystResults(const char *filein = "Chi_c_Syst_vDissertation", int nFiles = 6, const char *fileout = "Chi_c_Syst_Total_vDissertation.root", bool individualMerging = true) {
void mergeSystResults(const char *filein = "Chi_c_Jpsi_Syst_vDissertation", int nFiles = 6, const char *fileout = "Chi_c_Jpsi_Syst_Total_vDissertation.root", bool individualMerging = true) {
//void mergeSystResults(const char *filein = "Chi_c_Syst", int nFiles = 13, const char *fileout = "Chi_c_Syst_Total.root", bool individualMerging = true) {

	if (individualMerging == true) { // merge the individual files for each set as defined in ChiFitterInit.C
		for (int i = 0; i < nFiles; i++)
		{
			TString fileNameMerged = TString::Format("%s_%i.root", filein, i);
			TString fileNameSearched = TString::Format("%s_%i_*.root", filein, i);
			string cmd = string("hadd ") + string(fileNameMerged) + string( " `ls ") + string(fileNameSearched) + string("`");
			cout << cmd << endl;
			system(cmd.c_str());
		}
	}


   TFile *fout = new TFile(fileout,"RECREATE");
   //   //fout->cd();
   //TDirectory *TreeDir = gDirectory->mkdir(TreeName);
   //TreeDir->cd();
   

   
   for (int i = 0; i < (nFiles); i++)
   {
	   TString fileNameIn = TString::Format("%s_%i.root",filein, i);
	   cout << fileNameIn << endl;
	   CopyFile(fileNameIn);
   }

   fout->Close();

}
