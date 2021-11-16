#include <function.h>
#include <beam.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main( int argc, char **argv )
{

    //prepare for rootfile
    double ini_x, ini_y, ini_z, ini_energy;//Initial incident particle information
    double reac_x, reac_y, reac_z, reac_energy;// The reaction point information
    double cm_ang;//The random emit angle in the c.m. system
    double part_energy, part_det_energy;//the energy of 7Li left the target surface, and the detected 7Li with Si detector
    double det_x, det_y;// the 7Li position in Detector surface
    double phi;
    int total_event;

    //Save project information in TREE
    TTree *tree = new TTree("tree", "elastic scattering");

    //Initial particle information
    tree->Branch("ini_x", &ini_x, "ini_x/D");
    tree->Branch("ini_y", &ini_y, "ini_y/D");
    tree->Branch("ini_z", &ini_z, "ini_z/D");
    tree->Branch("ini_energy", &ini_energy, "ini_energy/D");

    //The reaction point information:position, reaction energy, 
    tree->Branch("reac_x", &reac_x, "reac_x/D");
    tree->Branch("reac_y", &reac_y, "reac_y/D");
    tree->Branch("reac_z", &reac_z, "reac_z/D");
    tree->Branch("reac_energy", &reac_energy, "reac_energy/D");
    
	//The emmited particle of 7Li:emited angle in the c.m. system , Yield of the d(6He,7Li)n
    tree->Branch("cm_ang", &cm_ang, "cm_ang/D");

	//The emitted energy from target surface and the detected energy within the si detector
    tree->Branch("part_energy", &part_energy, "part_energy/D");
    tree->Branch("part_det_energy", &part_det_energy, "part_det_energy/D");

    //The 7Li particle position in detector surface
    tree->Branch("det_x", &det_x, "det_x/D");
    tree->Branch("det_y", &det_y, "det_y/D");


    //Set primary beam and read parameters
    Beam *beam_test = new Beam();
    beam_test->read_parameters();
    beam_test->print_cond();

    //position[0]: location_x (cm)
    //position[1]: location_y (cm)
    //position[2]: location_z (cm)
    //position[3: Flightlength
    //Inci_energy[0]:the initial energy
    double position[3],direction_Start[3],Emit_particle[4],position_Emit[3], det_position[3];
    double Inci_Energy[1];

    TH1F *h_strip_7Li = new TH1F( "h_strip_7LI", "", 100, 0., 10 );
    TH1F *h_strip_cm_ang = new TH1F( "h_strip_cm_ang", "", 1000, 0., 180 );

    //for each theta angle, Calculating the 7Li numbers which can be reach in detector surface
    for(int i =1; i<181; i++)
    {
        ini_x=0.0, ini_y=0.0, ini_z=0.0, ini_energy=0.0;
        reac_x=0.0, reac_y=0.0, reac_z=0.0, reac_energy=0.0;
        cm_ang=0.0;
        part_energy=0.0, part_det_energy=0.0;
        det_x=0.0, det_y=0.0;    

        total_event=beam_test->generate_event(i);
        cout<<"angle in the c.m. system= "<<i<<endl;
        for(int event = 1; event < total_event; event++)
        {
            beam_test->generate_beam(position,direction_Start,Inci_Energy);
            ini_x = position[0];
            ini_y = position[1];
            ini_z = position[2];
            ini_energy = Inci_Energy[0];

            beam_test->reaction_loc_target(position,Inci_Energy);
		    reac_x = position[0];
		    reac_y = position[1];
		    reac_z = position[2];
		    reac_energy = Inci_Energy[0];     

            //  the angle phi was randomly filled among all the angles
            phi = 2 * M_PI * generate_standard();

		    beam_test->NuclearReaction(i, phi, direction_Start, Inci_Energy, Emit_particle);
            beam_test->leave_target(position, Emit_particle, position_Emit);
            part_energy = Emit_particle[3]/7.0;
            //cout<<part_energy<<endl;
            int flag;
            flag = beam_test->judge_detector(Emit_particle, position_Emit, det_position);
            if(flag == 1)
            {
                det_x = det_position[0];
                det_y = det_position[1];
                cm_ang = i;
                part_det_energy = beam_test->energy_detector(Emit_particle[3]); 
                part_det_energy = part_det_energy/7.0;
                h_strip_7Li->Fill(part_det_energy);
                h_strip_cm_ang->Fill(cm_ang);
                tree->Fill(); 
            }  
        }           
    }

    TString ofn = "../simulation.root";
    TFile *fout = new TFile(ofn, "recreate");
    tree->Write();
    h_strip_7Li->Write();
    h_strip_cm_ang->Write();
    fout->Close();

    cout << endl;
    cout << endl;
    cout << "<Created> ./simulation.root" << endl;
    cout << "...simulation completed!" << endl;
    return 0;
}
