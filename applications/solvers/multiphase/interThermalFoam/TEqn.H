surfaceScalarField kappaf=mixture.kappaf();

// Add Conv
volScalarField NoInterface (mixture.nearInterface()) ;
forAll( mesh.C (),celli)
  {
  if ( NoInterface [ celli ] > 0)
    NoInterface [ celli ] = 0;
  else
    NoInterface [ celli ] = 1;
  }
// End add Conv


fvScalarMatrix TEqn
(
fvm::ddt(rhoCp,T)
+fvm::div(rhoCpPhi,T)
-fvm::laplacian(kappaf,T)
==
// Radiative term
mixture.calcSourceRadiation(U ,T , RadiativeCoeff)
+ mixture.calcSourceForcedConvection(U , T, ConvectiveCoeff ) // Add the convective term
);
TEqn.relax();
TEqn.solve();
