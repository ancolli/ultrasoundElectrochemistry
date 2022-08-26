{   

    Info<< "min/max ap" << min(ap) << ' ' << max(ap) << endl;

    Ua = fvc::grad(ap)/(omega_f);//rho
    phiUa = linearInterpolate(Ua) & mesh.Sf();

    forAll(mesh.C(), celli)
    {
        if (ap[celli] > 103.0)//(ap[celli] > P_Blake.value()){ 
        {
        //bubbleRadialDynamics from https://www.sciencedirect.com/science/article/pii/S1350417722001535
            k_r_sq[celli] = sqr(omega_f.value()/c.value())+A.value();
            k_i_sq[celli] = -B.value()*(172.1*Foam::exp(0.036*ap[celli])/(172.1*195.4+Foam::exp(0.036*ap[celli]))            
            +22.0*Foam::exp(0.02*ap[celli])/(22.0*10.0+Foam::exp(0.02*ap[celli])))/Foam::pow(ap[celli],2.0);
        }

        else {
            k_r_sq[celli] = sqr(omega_f.value()/c.value());
            k_i_sq[celli] = 0.0;//
        }

    }
    //volScalarField ap_old = max(ap.oldTime(), dimensionedScalar("min_ap", dimensionSet(1,-1,-2,0,0,0,0), 1.0));
    //k_r_sq == sqr(omega/c)-A*((2/ap_old)-1/pow(ap_old,2)*ap);
    //k_i_sq == -B*((2/ap_old)-1/pow(ap_old,2)*ap);
    
    Info<< "min/max k_r_sq_" << min(k_r_sq) << ' ' << max(k_r_sq) << endl;
    Info<< "min/max k_i_sq_" << min(k_i_sq) << ' ' << max(k_i_sq) << endl;
      
    p_Re.storePrevIter();
    p_Im.storePrevIter();
   
    fvScalarMatrix ReEqn
    (
        fvm::laplacian(p_Re)
        + fvm::SuSp(k_r_sq,p_Re)
        == k_i_sq*p_Im
        //+ fvm::SuSp(k_r_sq-k_i_sq*p_Im.oldTime()/p_Re.oldTime(),p_Re)//+ fvm::SuSp(k_P, p_Re)//
        
    );
    
    fvScalarMatrix ImEqn
    (
        fvm::laplacian(p_Im)     
        + fvm::SuSp(k_r_sq,p_Im)
        == - k_i_sq*p_Re//Bi
        //+ fvm::SuSp(k_r_sq+k_i_sq*p_Re.oldTime()/p_Im.oldTime(),p_Im)//+ fvm::SuSp(k_P, p_Im)//
        
    );
        
    p_Re.relax();
    p_Im.relax();
    
    ReEqn.solve();
    ImEqn.solve();
   
    p_Re.correctBoundaryConditions();
    p_Im.correctBoundaryConditions();
    
    ap = (sqrt((pow(p_Re, 2) + pow(p_Im, 2))/2));
    
    ap.correctBoundaryConditions();
}



    
 

    
 



