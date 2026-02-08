// Simple example burner model for use with new NPSS installations
// SET THERMODYNAMIC PACKAGE
setThermoPackage("Janaf");
// SPECIFY MODEL ELEMENTS
// Start the flow of air
Element FlowStart FsAir {
 Pt {value=14.7; units="psia";} // psi
 Tt = 519.; // R
 W = 10.; // lbm/s
 WAR = 0; // Water to Air ratio
}
// Start the flow of fuel
Element FuelStart Fus {
 Wfuel = 0.3; // lbm/s, fuel flow rate (Used ONLY when Burner switchBurn = WFUEL)
 LHV = 18000; // BTU/lbm, user input fuel lower heating value (LHV). Default is 18400 BTU/lbm
}
// Burn the fuel and air
Element Burner Brn {
 dPqP_dmd = 0.05; // user input friction relative pressure drop (Pin - Pout)/Pin
 effBase = 0.98; // user input burner adiabatic efficiency
 // The value for switchBurn determines how burner fuel flow rate is set
 switchBurn = "WFUEL"; // WFUEL, use Fu_I.Wfuel
}
// End the flow of air
Element FlowEnd FeAir;
// SET ELEMENT LINKAGES
linkPorts( "FsAir.Fl_O" , "Brn.Fl_I" , "F030" );
linkPorts( "Brn.Fl_O" , "FeAir.Fl_I" , "F040" );
linkPorts( "Fus.Fu_O" , "Brn.Fu_I" , "Fu010" );
// SET SOLUTION PARAMETERS
// In this example, let NPSS determine how much fuel flow is necessary
// to achieve a specific burner temperature
// Define fuel flow as independent variable
Independent ind_Wfuel {
 varName = "Fus.Wfuel";
 indepRef = "0.5";
 dxLimit = 0.2;
 dxLimitType = "FRACTIONAL";
 perturbation = 0.01;
 perturbationType = "FRACTIONAL";
 description = "vary the fuel flow to achieve the desired burner temp";
}
// Define burner temperature as dependent variable
Dependent dep_F040Tt {
 eq_rhs = "Tt4_in";
 eq_lhs = "F040.Tt";
 eq_Ref = "Tt4_in";
 description = "desired burner temperature";
}
// Add independent/dependent pair to solver
solver.addDependent("dep_F040Tt");
solver.addIndependent("ind_Wfuel");
// RUN SOLUTION
// Declare Tt4_in and set desired value of burner temperature
real Tt4_in { value = 2500.; units = "R"; };
// Execute solution
run();
// DISPLAY DESIRED RESULTS TO SCREEN
cout << "F040.Tt = " << F040.Tt << " " << F040.Tt.units << endl;
cout << "Wfuel = " << Brn.Wfuel << " " << Brn.Wfuel.units << endl;