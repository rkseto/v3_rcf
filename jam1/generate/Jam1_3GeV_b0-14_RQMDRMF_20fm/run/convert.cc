#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TTree.h"
#include "TFile.h"
#include "TJam.h"

using TJAM::TJam;
const int split_step = 1000; //How to split events into root files.

int EventLoop(TJam&);

// The callback function for time evolution recording
void record_time(TJam& jam, float atime);

int timesel[] = {5, 10, 15, 20, 30}; // Record at this time (fm/c)
int lentime = sizeof(timesel) / sizeof(int);

const int maxv = 30000;
struct RecordField{
    int Npart, mul;
    float b, atime;
    int pid[maxv];
    int ks[maxv];
    float rx[maxv], ry[maxv], rz[maxv], rt[maxv];
    float frx[maxv], fry[maxv], frz[maxv], frt[maxv];
    float px[maxv], py[maxv], pz[maxv];
    //float E[maxv], mass[maxv], charge[maxv];
    //float pt[maxv], ycm[maxv], eta[maxv];
} RecField;

TFile *TimeRecord[100];
TTree *TimeTree[100];

int main(void) {
    using namespace std;
    //Generate random number
    int rd;
    ifstream fin("rdlist");
    fin>>rd;
    fin.close();

    //Ready to emit jam
    TJam jam;
    jam.mstc(1)=rd;     //Random number seed
    jam.mevent=200;    //Event amount
    //jam.mevent=1000; //Event amount
    jam.bmin=0;       //Minimum impact parameter b
    jam.bmax=-13.42;     //Maximum impact parameter b
    jam.cwin="3gev";  //Incident energy, Ecm here. See /src/jam/README
    jam.frame="collider";   //Compution frame, cm, nn, lab or collider
    jam.proj="197Au";   //Project
    jam.targ="197Au";   //Target
    //jam.dt=40;
    //jam.nstep=1;
    jam.dt=0.2;       //Time step in fm/c
    jam.nstep=100;
    jam.mstc(6)=204;  //Mean Field, RQMD/S
    //jam.mstc(42)=0;     //Allow weak decy
    //jam.mstc(51)=0;       //Only baryon-baryon collisions
    //jam.mstc(34)=1;       //No spectators interaction
    //jam.mstc(17)=1;       //only inelastic collisions are generated
    //jam.record_time_evolution = true;
    //jam.record_time_function = record_time;
	//jam.mstc(59)=22;	//attractive orbit
    jam.quiet_file_out(3); // 1-4 are: '0','JAMRUN.DAT','JAMINFO.DAT','JAMMULTI.DAT'
    lentime = -1;
    //Emit jam
    jam.jaminit();

    //Splited Event Loop, to avoid memory leaks
    for(int iev = 1; iev <= jam.mevent; ) {
        iev = EventLoop(jam);
    }//End Event Loop

    //Finish simulation, output
    jam.jamfin();
    return 0;
}

int EventLoop(TJam& jam){
    using namespace std;
    //<step> events will be recorded in one .root file;
    const int step = split_step;
    int start = jam.iev() + 1;
    int end = start + step - 1;
    end = end > jam.mevent ? jam.mevent : end;
    static char RootName[255], TreeName[255];
    int i, nv, mul, iev;
    //New Tree
    sprintf(RootName, "event_%d.root", end);
    TTree EventTree("jam", "jam");
    //Init Branch
    EventTree.Branch("mul", &RecField.mul, "mul/I");
    EventTree.Branch("b", &RecField.b, "b/F");
    EventTree.Branch("Npart", &RecField.Npart, "Npart/I");
    EventTree.Branch("ks", &(RecField.ks), "ks[mul]/I");
    EventTree.Branch("pid", RecField.pid, "pid[mul]/I");
    EventTree.Branch("x", RecField.rx, "x[mul]/F");
    EventTree.Branch("y", RecField.ry, "y[mul]/F");
    EventTree.Branch("z", RecField.rz, "z[mul]/F");
    EventTree.Branch("t", RecField.rt, "t[mul]/F");
    EventTree.Branch("frx", RecField.frx, "frx[mul]/F");
    EventTree.Branch("fry", RecField.fry, "fry[mul]/F");
    EventTree.Branch("frz", RecField.frz, "frz[mul]/F");
    EventTree.Branch("frt", RecField.frt, "frt[mul]/F");
    EventTree.Branch("px", RecField.px, "px[mul]/F");
    EventTree.Branch("py", RecField.py, "py[mul]/F");
    EventTree.Branch("pz", RecField.pz, "pz[mul]/F");
    //EventTree.Branch("E", RecField.E, "E[mul]/F");
    //EventTree.Branch("mass", RecField.mass, "mass[mul]/F");
    //EventTree.Branch("charge", RecField.charge, "charge[mul]/F");
    //EventTree.Branch("pt", RecField.pt, "pt[mul]/F");
    //EventTree.Branch("ycm", RecField.ycm, "ycm[mul]/F");
    //EventTree.Branch("eta", RecField.eta, "eta[mul]/F");

    for(iev = start; iev <= end; ++iev) {
        if(iev % step == 0) cout << "# Event = " << iev << endl;
        //Generate an event
        jam.jamevt(iev);

        //nv = Total Particle Number in this event
        RecField.mul = jam.mul();
        RecField.b = jam.pard(2);
        //Particle Loop
        jam.write_cldata(&RecField.Npart, RecField.ks, RecField.pid,
                RecField.rx, RecField.ry, RecField.rz, RecField.rt,
                RecField.frx, RecField.fry, RecField.frz, RecField.frt,
                RecField.px, RecField.py, RecField.pz);

        //if (RecField.Npart == 0) { // skip empty
        //    --iev;
        //    continue;
        //}
        EventTree.Fill();
    }//End Even Loop
    //Write to .root
    TFile *EventRecord = new TFile(RootName, "RECREATE");
    EventRecord->cd();
    EventTree.Write();
    //EventRecord.Write();
    EventRecord->Close();
    delete EventRecord;

    //Return where the next event loop start
    return iev;
}

void record_time(TJam &j, float atime) {
    // This function is call by TJam to record the time evolution
    // Parameters tell it who's calling and when(TJam&, float)
    using namespace std;
    bool drop = true;
    int i;
    for (i = 0; i < lentime; ++i)
        if (int(atime) == timesel[i]) {
            drop = false;
            break;
        }
    if (drop) return;
    cout << "aTime = " << atime << " timestep = ";
    cout << j.mstd(23) << endl;
    RecField.mul = j.mul();
    // useless, because mul will be changed in write_time(...)
    // since nv is not the total number of particle to record
    // subroutine will jump particles dead or not formed

    /*
    j.write_time(&(RecField.atime), &(RecField.Npart), &(RecField.mul),
            RecField.ks, RecField.pid, RecField.rx, RecField.ry, RecField.rz,
            RecField.rt, RecField.rft,
            RecField.px, RecField.py, RecField.pz,
            RecField.E, RecField.mass, RecField.charge,
            RecField.pt, RecField.ycm, RecField.eta);
    */

    RecField.b = j.pard(2);

    TimeTree[i]->Fill();
}
