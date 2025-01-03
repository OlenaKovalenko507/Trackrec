#include <TSystem.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TText.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>

using namespace std;

// Declare the hit structure before any function that uses it
struct hit {
    int x_i;  // x-coordinate of the hit
    int y_i;  // y-coordinate of the hit
    double r_i;  // radius of the hit
};

struct track {
    double x; 
    double y; 
    double r; 
}; 

const int n1 = 10;
const int n2 = 10;
const int n_points = 10000; // number of points used to approximate the unit circle 
double delta = 0.5; // size of the bin 
const double zoom_scale = 10; // zoom-in scale between iterations 
int num_of_iterations = 4; 
double x_min = -10; 
double y_min = -10; 
double x_max = 10 + n1; 
double y_max = 10 + n2; 
double z_min = 1; 
double z_max = 10;

void CircularEventGenerator(double tracker_cells[][n2], double x0, double y0, double r0);
track CreateConeCrossections(vector<hit>& event, double z_max, int n_points);  // Pass by reference
void CircularEventGenerator( double tracker_cells[][n2], double x0, double y0, double r0, double phi_min, double phi_max);

int CircTrack() {
    //////////////////////////
    // setting up tracker hits
    //////////////////////////

    double tracker_cells[n1][n2] = {0.0};

    // define the event
    CircularEventGenerator(tracker_cells, 6.232, 5.23, 3.6);

    // create a vector containing all events
    vector<hit> event;
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            if (tracker_cells[i][j] != 0.0) {
                hit h = {i, j, tracker_cells[i][j]};
                event.push_back(h);  // push the hit into the vector
            }
        }
    }

    // total number of hits
    int total_no_hits = event.size();
    cout << "number of hits: " << total_no_hits << endl;



    // Call CreateConeCrossections
    track final_track = CreateConeCrossections(event, z_max, n_points);  

    TEllipse* big_circle = new TEllipse(final_track.x, final_track.y, final_track.r, final_track.r); 
    
    // creating elipses for vizualization purposes
	TEllipse *el[total_no_hits];
	for(int i = 0; i < total_no_hits; i++)
	{
		el[i] = new TEllipse(event[i].x_i, event[i].y_i, event[i].r_i, event[i].r_i);	
	}
    TCanvas *can2 = new TCanvas("can1", "canvas", 1000, 1000);
    big_circle->Draw("SAME"); 
	can2->Range(-2, -2, n1+2, n2+2);
	for(int i = 0; i < total_no_hits; i++)
	{
		el[i]->SetLineWidth(2);
   		el[i]->Draw("SAME");
	}
    
	can2->SaveAs("circular_track.png");

    return 0;
}

// function that generates hits that lie on a circular trajectory
void CircularEventGenerator(double tracker_cells[][n2], double x0, double y0, double r0) {
    double distance;
    double deltaxy;

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            deltaxy = (i - x0) * (i - x0) + (j - y0) * (j - y0);
            distance = abs(r0 - sqrt(deltaxy));
            if (tracker_cells[i][j] == 0.0 && distance < 0.5)
                tracker_cells[i][j] = distance;
        }
    }
}

void CreateCircle(double x0, double y0, double r, int n_points, TH2F* tracker_hits) {
    int points_on_circle = int(n_points * r);  // Calculate the number of points for this circle based on its radius
    for (int i = 0; i < points_on_circle; ++i) {
        double theta = 2.0 * M_PI * i / points_on_circle;  // Angle from 0 to 2*pi
        double x = r * cos(theta) + x0;
        double y = r * sin(theta) + y0;

        // Check which bin (cell) the point (x, y) falls into and fill the histogram
        int bin_x = tracker_hits->GetXaxis()->FindBin(x);
        int bin_y = tracker_hits->GetYaxis()->FindBin(y);
        
        tracker_hits->Fill(x, y);
        
    }
}

// Create cross sections of the cones for each tracker hit
track CreateConeCrossections(vector<hit>& event, double z_max, int n_points) {
    double x; 
    double y; 
    double r; 
    int x_peak; 
    int y_peak; 
    int z_peak; 
    TH2F* cross_sec; 

    //size of 1 bin in the histogram (changes with each iteration)
    for (int iteration = 0; iteration < num_of_iterations; iteration++) {
    // positions of the maximum (changes with each iteration)
        int global_max_bin = 0; //Global maximum bin
        double global_max_value = 0; //Maximum values across all diagrams 

        double max_z_value; 

        cout << "x_max: " << x_max << endl; 
        cout << "x_min: " << x_min << endl; 
        cout << "y_min: " << y_min << endl; 
        cout << "y_max: " << y_max << endl; 

    // Create the 2D histogram (tracker_hits) for this cross-section
        for (double z = z_min; z <= z_max; z += delta) { 
            cross_sec = new TH2F(Form("tracker_z_/iteration%f %d", z, iteration), Form("Cross Sections at z = %f %d; x; y", z, iteration), (x_max - x_min)/delta, x_min, x_max, (y_max - y_min)/delta, y_min, y_max); // ADD ITERATION NUMBER TO THE NAME 
            
            for (auto& itr : event) {  // Reference to avoid copying
                CreateCircle(itr.x_i, itr.y_i, z + itr.r_i, n_points, cross_sec);  // Pass all required arguments
                CreateCircle(itr.x_i, itr.y_i, z - itr.r_i, n_points, cross_sec);  // Pass all required arguments
            }

            // Create a canvas to draw the histogram
            TCanvas* c = new TCanvas(Form("c_z_%f %d", z, iteration), Form("Cross Sections at z = %f %d", z, iteration), 800, 600);
            cross_sec->Draw("COLZ");

            // Loop over the bins and display the number of hits
            for (int ix = 1; ix <= cross_sec->GetNbinsX(); ++ix) {
                for (int jy = 1; jy <= cross_sec->GetNbinsY(); ++jy) {
                    int bin_content = cross_sec->GetBinContent(ix, jy);
                    if (bin_content > 0) {
                        double x = cross_sec->GetXaxis()->GetBinCenter(ix);
                        double y = cross_sec->GetYaxis()->GetBinCenter(jy);
                        TString label = Form("%d", bin_content);
                        TText *text = new TText(x, y, label);
                        text->SetTextAlign(22);  // Center alignment
                        text->SetTextSize(0.025);
                        //text->Draw();
                    }
                }
            }
            c->Update();  // Update the canvas to show all drawings
            // Find the bin with the maximum content for this histogram  
            int local_max_bin = cross_sec->GetMaximumBin();
            double local_max_value = cross_sec->GetBinContent(local_max_bin); 

            // Compare it with the global maximum and update if necessary
            if (local_max_value > global_max_value) {
                global_max_value = local_max_value;
                global_max_bin = local_max_bin;
                max_z_value = z;
                r = z; 
            }
        

        }
    //find values of x and y for the global maximum  
	cross_sec->GetBinXYZ(global_max_bin, x_peak, y_peak, z_peak);
    x = cross_sec->GetXaxis()->GetBinCenter(x_peak);
    y = cross_sec->GetYaxis()->GetBinCenter(y_peak);


    cout << "Global Maximum Value: " << global_max_value << endl;
    cout << "Global Maximum Bin: " << global_max_bin << endl;
    cout << "At z value: " << max_z_value << endl;

    double x_length = x_max - x_min; 
    double y_length = y_max - y_min; 
    double z_length = z_max - z_min; 
    x_min = x - 0.5*x_length/zoom_scale; 
    x_max = x + 0.5*x_length/zoom_scale; 
    y_min = y - 0.5*y_length/zoom_scale; 
    y_max = y + 0.5*y_length/zoom_scale; 
    z_min = max_z_value - 0.5*z_length/zoom_scale; 
    z_max = max_z_value + 0.5*z_length/zoom_scale; 
    delta = delta/zoom_scale;
    }
    track big_circle = {x, y, r}; 
    return big_circle; 
}

void CircularEventGenerator( double tracker_cells[][n2], double x0, double y0, double r0, double phi_min, double phi_max)
{
	double distance;
	double deltaxy;
	double phi;
	
	for(int i = 0; i < n1; i++)
		for(int j = 0; j < n2; j++)
		{
			deltaxy =  ((double)i-x0)*((double)i-x0) + ((double)j-y0)*((double)j-y0);
			phi = 180*(j-y0)*acos( (i-x0)/sqrt(deltaxy) ) / (M_PI * abs(j-y0));	
			if(phi > phi_min && phi < phi_max)
			{
				distance = abs(r0 - sqrt(deltaxy));
				if( tracker_cells[i][j] == 0.0 && distance < 0.5 )
					tracker_cells[i][j] = distance;
			}
		}
}
