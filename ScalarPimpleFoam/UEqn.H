 
    // Solve the momentum equation

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevSigma(U)
     ==
        fvModels.source(U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvConstraints.constrain(U);
    }
