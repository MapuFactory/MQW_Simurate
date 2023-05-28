classdef Constant
    properties
        ML
        MSTAR
        ELEC
        DIEELECSTAR
        HBAR
        HH
        KB
        TEMP
        KT
        II
        VI
        DELTA
        DELTAE
        DX
        OUTPUT
        NI
        EG_SI
        BASE
        pBASE
        DIE_CAF2
        EAFF_CAF2
        BAR_CAF2
        MASS_CAF2
        DIE_CDF2
        EAFF_CDF2
        BAR_CDF2
        MASS_CDF2
        DIE_iSI
        EAFF_iSI
        BAR_iSI
        MASS_iSI
        NUMVALLY_SI
        MASS_SI_Z
        MASS_SI_XY
        DIE_AL
        EAFF_AL
        BAR_AL
        MASS_AL
        EF_AL
        DIE_AU
        EAFF_AU
        BAR_AU
        MASS_AU
        EF_AU
        DIE_pSI
        EAFF_pSI
        BAR_pSI
        MASS_pSI
        DIE_pAL
        EAFF_pAL
        BAR_pAL
        MASS_pAL
        EF_pAL
        DIE_pCAF2
        EAFF_pCAF2
        BAR_pCAF2
        MASS_pCAF2
        DIE_SIO2
        EAFF_SIO2
        BAR_SIO2
        MASS_SIO2
        NUMVALLY_SI_Z
        NUMVALLY_SI_XY
        LBAR
        WELL
        RBAR
        DIM
        MAX_LAYER
        DIV
        EMAX
        AAA
        RTD
        ALL
        dx
    end
    
    methods
        function obj = Constant()
            obj.ML = 0.31e-9;
            obj.MSTAR = 9.109e-31;
            obj.ELEC = 1.602e-19;
            obj.DIEELECSTAR = 8.854e-12;
            obj.HBAR = 1.054e-34;
            obj.HH = 6.626e-34;
            obj.KB = 1.38e-23;
            obj.TEMP = 3.0e+2;
            obj.KT = obj.KB * obj.TEMP;
            obj.II = obj.ELEC * obj.ELEC * obj.MSTAR * obj.KB / 2 / (pi*pi) / obj.HBAR / obj.HBAR / obj.HBAR;
            obj.VI = 0.6;
            obj.DELTA = 1e-5;
            obj.DELTAE = 1e-9;
            obj.DX = 10;
            obj.OUTPUT = 1;
            obj.NI = 1.45e16;
            obj.EG_SI = 1.12;
            obj.BASE = 4.05;
            obj.pBASE = -5.17;
            obj.DIE_CAF2 = 6.76;
            obj.EAFF_CAF2 = 4.05-1.0;
            obj.BAR_CAF2  = -(obj.EAFF_CAF2-obj.BASE);
            obj.MASS_CAF2 = 1;
            obj.DIE_CDF2 = 8.83;
            obj.EAFF_CDF2 = (4.05+0.6);
            obj.BAR_CDF2 = -(obj.EAFF_CDF2-obj.BASE);
            obj.MASS_CDF2 = 0.4;
            obj.DIE_iSI = (11.8);
            obj.EAFF_iSI = 4.05;
            obj.BAR_iSI = -(obj.EAFF_iSI-obj.BASE);
            obj.MASS_iSI = 0.26;
            obj.NUMVALLY_SI = 6;
            obj.MASS_SI_Z = 0.98;
            obj.MASS_SI_XY = 0.19;
            obj.DIE_AL = 0;
            obj.EAFF_AL = 15.83;
            obj.BAR_AL = -(obj.EAFF_AL-obj.BASE);
            obj.MASS_AL = 1;
            obj.EF_AL = (11.63);
            obj.DIE_AU = 0;
            obj.EAFF_AU = 10.33;
            obj.BAR_AU = -(obj.EAFF_AU-obj.BASE);
            obj.MASS_AU = 1;
            obj.EF_AU = 5.51;
            obj.DIE_pSI = (11.8);
            obj.EAFF_pSI = -(4.05+1.12);
            obj.BAR_pSI = -(obj.EAFF_pSI-obj.pBASE);
            obj.MASS_pSI = 0.55;
            obj.DIE_pAL = 0;
            obj.EAFF_pAL = 0;
            obj.BAR_pAL = -(obj.EAFF_pAL-obj.pBASE);
            obj.MASS_pAL = 1;
            obj.EF_pAL = 4.37;
            obj.DIE_pCAF2 = 6.76;
            obj.EAFF_pCAF2 = -(4.05+1.12+1);
            obj.BAR_pCAF2 = -(obj.EAFF_pCAF2-obj.pBASE);
            obj.MASS_pCAF2 = 1;
            obj.DIE_SIO2 = 3.9;
            obj.EAFF_SIO2 = (4.05-2.9);
            obj.BAR_SIO2 = -(obj.EAFF_SIO2-obj.BASE);
            obj.MASS_SIO2 = 0.26;
            obj.NUMVALLY_SI_Z = 2;
            obj.NUMVALLY_SI_XY = 4;
            obj.LBAR = 2;
            obj.WELL = (obj.LBAR+1);
            obj.RBAR = (obj.LBAR+2);
            obj.DIM = 2;
            obj.MAX_LAYER = 100;
            obj.DIV = 10;
            obj.EMAX = 10;
            obj.AAA = 0;
            obj.RTD = 2;
            obj.ALL = 1;
            obj.dx=obj.ML/obj.DX;	%/* DXは分割数	MLはこのプログラムのz軸の基本単位。そのMLをさらにDX分割する。*/
        end
    end
end