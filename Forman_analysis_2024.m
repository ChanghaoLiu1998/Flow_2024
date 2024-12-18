% -------------------------------------------------------------------------
% preliminaries
% -------------------------------------------------------------------------
clc;
close all;
clear;


    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    filename = sprintf('mytess1.tess');
    disp(['>> Processing ', filename]);
    
    % ---------------------------------------------------------------------
    % create Voronoi complex structure
    % ---------------------------------------------------------------------
    vcmplx = setup_complex(4);
    vcmplx(1).name='vertex'; vcmplx(2).name='edge';
    vcmplx(3).name='face';   vcmplx(4).name='polyhedron';
    
    % ---------------------------------------------------------------------
    % read Neper .tess file for Voronoi complex
    % ---------------------------------------------------------------------
    vf = 0.214;
    disp(['>> read Neper .tess file (', filename, ') for Voronoi complex']);
    [vcmplx] = read_tess(filename, vcmplx);
    
    % ---------------------------------------------------------------------
    % get Voronoi geometry
    % ---------------------------------------------------------------------
    disp('>> get Voronoi geometry');
    vcmplx = vgeo(vcmplx);

    % ---------------------------------------------------------------------
    % build Voronoi connectivity
    % ---------------------------------------------------------------------
    disp('>> build Voronoi connectivity');
    vcmplx = vcon(vcmplx);
    
    % ---------------------------------------------------------------------
    % create Forman complex structure
    % ---------------------------------------------------------------------
    disp('>> create Forman complex structure');
    fcmplx = setup_complex(4);
    fcmplx(1).name='vertex'; fcmplx(2).name='edge';
    fcmplx(3).name='face';   fcmplx(4).name='polyhedron';
    
    % ---------------------------------------------------------------------
    % setup Forman complex from Voronoi complex
    % ---------------------------------------------------------------------
    disp('>> setup Forman complex from Voronoi complex');
    fcmplx = v2f(fcmplx, vcmplx);
    
    % ---------------------------------------------------------------------
    % build Forman connectivity - boundary operator fBO_1
    % ---------------------------------------------------------------------
    disp('>> build Forman connectivity');
    fcmplx = fcon(fcmplx);
    fBO_1 = build_BO(fcmplx, 1);
    
    % ---------------------------------------------------------------------
    % get Forman geometry
    % ---------------------------------------------------------------------
    disp('>> get Forman geometry');
    fcmplx = fgeo(fcmplx);
    
    % ---------------------------------------------------------------------
    % get Forman coboundary star - fCBS_1
    % ---------------------------------------------------------------------
    disp('>> get Forman coboundary star');
    ifs1 = spdiags(1./fcmplx(1).fvol, 0, fcmplx(1).num(1).val, fcmplx(1).num(1).val);
    fs2  = spdiags(fcmplx(2).fvol, 0, fcmplx(2).num(2).val, fcmplx(2).num(2).val);
    fCBS_1 = -ifs1 * fBO_1' * fs2;
    
    % ---------------------------------------------------------------------
    % musical iso (flat, sharp)
    % ---------------------------------------------------------------------
    fcmplx = musical_iso_forman(fcmplx, fBO_1);
    
    % ---------------------------------------------------------------------
    % set up graphite model for diffusion coefficients
    % ---------------------------------------------------------------------
    model;
    
    % ---------------------------------------------------------------------
    % solve pre-oxidised graphite
    % ---------------------------------------------------------------------
    solve_forman;
    
  
   