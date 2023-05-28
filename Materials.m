classdef Materials < handle
	properties
		die;		%/* 誘電率					*/
		bar;		%/* バンド不連続				*/
		mass;	%/* 有効質量					*/
		massxy;	%/* 横方向の有効質量			*/
		ef;		%/* フェルミエネルギー		*/
		valley;		%/* 伝導帯の谷の数			*/
		name;		%/* 物質の名前				*/
		cond;
		Ef;		%フェルミエネルギーの設定。金属は物性値、半導体は教科書の式から計算、絶縁物には対応していない。
		base;
		Q;
		Work;
		divnum;
		NX;
		smt;
        d;
	end

	methods

		function obj = Materials(materialName, ML, Q, NX)
			Ed = -0.05;
			base = 0;		%/* 伝導帯エネルギーの基準。大抵はSiの電子親和力	*/
            const = Constant();
			if nargin == 0
				obj.die = 0;
				obj.bar = 0;
				obj.mass = 0;
				obj.massxy = 0;
				obj.ef = 0;
				obj.valley = 0;
				obj.name = 0;
				obj.cond = 0;
				obj.Ef = 0;
				obj.Q = 0;
				obj.base = const.pBASE;
				obj.Work = -obj.Ef - obj.bar + obj.base;
				obj.divnum = 0;
				obj.NX = 0;
				obj.smt = 0;
                obj.d = 0;
			else
				if strcmp(materialName, 'Al')
					obj.name = 'Al';
					obj.die = const.DIE_AL * const.DIEELECSTAR ;
					obj.bar = const.BAR_AL ;
					obj.mass = const.MASS_AL ;
					obj.massxy = const.MASS_AL;
					obj.valley = 1;
					obj.cond = 0;
					obj.Ef = -0.56+const.EF_AL;
					obj.Q = Q;
					obj.base = const.BASE; % //-EG_SI/2+KT*log(Q[i]*1e6/NI)/ELEC;base=BASE;break;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 0;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'n-Si.Sub.')
					obj.name = 'n-Si.Sub.';
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = const.MASS_iSI;
					obj.massxy = const.MASS_iSI;
					obj.valley = const.NUMVALLY_SI;
					obj.cond=1;
					obj.Q = Q;
					obj.Ef = -const.EG_SI/2+const.KB*const.TEMP*log((obj.Q*1e6)/const.NI)/const.ELEC;
					obj.base = const.BASE;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 1;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'n-Si')
					obj.name="n-Si";
					obj.die = const.DIE_iSI  *const. DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = const.MASS_iSI;
					obj.massxy = const.MASS_iSI;
					obj.valley = const.NUMVALLY_SI;
					obj.cond=1;
					obj.Q = Q;
					obj.Ef = -const.EG_SI/2+const.KB*const.TEMP*log((obj.Q*1e6)/const.NI)/const.ELEC;
					obj.base = const.BASE;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 2;
                    obj.d = ML*const.ML;
				% elseif materialName, 'p-Si.Sub.'
				% 	obj.name="p-Si.Sub.";
				% 	obj.die=DIE_iSI  * DIEELECSTAR;
				% 	obj.bar=BAR_iSI;
				% 	obj.mass=MASS_iSI;
				% 	obj.massxy=MASS_iSI;
				% 	obj.valley=NUMVALLY_SI;
				% 	obj.cond=2;
				%	obj.smt = 3;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'p-Si')
					obj.name = "p-Si";
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = const.MASS_iSI;
					obj.massxy = const.MASS_iSI;
					obj.valley = const.NUMVALLY_SI;
					obj.cond = 2;
					obj.Q = Q;
					nv = 2*pow(const.MSTAR*const.KB*const.TEMP/const.HBAR/const.HBAR/pi/2.0, 1.5 ) * (pow(0.16,1.5) + pow(0.49,1.5));
					obj.Ef = Ed + const.KB*const.TEMP*log(-1.0/4.0 + pow(1+8*obj.Q*1e6*exp(-Ed/const.KB/const.TEMP)/nv, 0.5)/4)/const.ELEC;
					obj.base = const.pBASE;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 4;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'CaF2')
					obj.name = "CaF2";
					obj.die = const.DIE_CAF2 * const.DIEELECSTAR;
					obj.bar = const.BAR_CAF2;
					obj.mass = const.MASS_CAF2;
					obj.massxy = const.MASS_CAF2;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 0;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 5;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'CdF2')
					obj.name = "CdF2";
					obj.die = const.DIE_CDF2 * const.DIEELECSTAR;
					obj.bar = const.BAR_CDF2;
					obj.mass = const.MASS_CDF2;
					obj.massxy = const.MASS_CDF2;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 0;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 6;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'i-Si')
					obj.name = "i-Si";
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = const.MASS_iSI;
					obj.massxy = const.MASS_iSI;
					obj.valley = const.NUMVALLY_SI;
					obj.cond = 3;
					obj.Ef = -0.56;
					obj.Q = Q;
					obj.base = const.BASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 7;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'Au')
					obj.name = "Au";
					obj.die = const.DIE_AU   * const.DIEELECSTAR;
					obj.bar = const.BAR_AU ;
					obj.mass = const.MASS_AU;
					obj.massxy = const.MASS_AU;
					obj.valley = 1;
					obj.cond = 0;
					obj.Ef = const.EF_Au;
					obj.Q = Q;
					obj.base = const.BASE;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 8;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'nSub.CaF2')
					obj.name = "nSub.CaF2";
					obj.die = const.DIE_CAF2 * const.DIEELECSTAR;
					obj.bar = const.BAR_CAF2;
					obj.mass = const.MASS_CAF2;
					obj.massxy = const.MASS_CAF2;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = -0.3;
					obj.Q = Q;
					obj.base = const.BASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 9;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'p-Si.Sub.')
					obj.name = "p-Si.Sub.";
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_pSI;
					obj.mass = MASS_pSI;
					obj.massxy = const.MASS_pSI;
					obj.valley = 1;
					obj.cond = 3;
					obj.Q = Q;
					nv = 2*pow(const.MSTAR*const.KB*TEMP/const.HBAR/const.HBAR/pi/2.0, 1.5 ) * (pow(0.16,1.5) + pow(0.49,1.5));
					obj.Ef = Ed + const.KB*const.TEMP*log(-1.0/4.0 + pow(1+8*obj.Q*1e6*exp(-Ed/const.KB/const.TEMP)/nv, 0.5)/4)/const.ELEC;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 10;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'pCaF2')
					obj.name = "pCaF2";
					obj.die = const.DIE_pCAF2* const.DIEELECSTAR;
					obj.bar = const.BAR_pCAF2;
					obj.mass = const.MASS_pCAF2;
					obj.massxy = const.MASS_pCAF2;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 0;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 11;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'pAl')
					obj.name = "pAl";
					obj.die = const.DIE_pAL  * const.DIEELECSTAR;
					obj.bar = const.BAR_pAL;
					obj.mass = const.MASS_pAL;
					obj.massxy = const.MASS_pAL;
					obj.valley = 1;
					obj.cond = 0;
					obj.Ef = const.EF_pAL;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 12;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, '5MLCaF2')
					obj.name = "5MLCaF2";
					obj.die = const.DIE_CAF2 * const.DIEELECSTAR;
					obj.bar = 1.7;
					obj.mass = 0.7;
					obj.massxy = 0.7;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 0;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 13;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, '仮想AL')
					obj.name = "仮想AL";
					obj.die = const.DIE_AL   * const.DIEELECSTAR;
					obj.bar = -5.17  ;
					obj.mass = 1;
					obj.massxy = 1;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 4.37;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 14;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, '仮想CaF2')
					obj.name = "仮想CaF2";
					obj.die = const.DIE_CAF2 * const.DIEELECSTAR;
					obj.bar = 7;
					obj.mass = 1;
					obj.massxy = 1;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 0;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 15;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, '仮想p-Si')
					obj.name = "仮想p-Si";
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = 0.55;
					obj.massxy = 0.55;
					obj.valley = 1;
					obj.cond = 3;
					obj.Q = Q;
					obj.Ef = -const.EG_SI/2 + const.KB*const.TEMP*log(log.Q*1e6/const.NI)/const.ELEC;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 16;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'SiO2')
					obj.name = "SiO2";
					obj.die = const.DIE_SIO2 * const.DIEELECSTAR;
					obj.bar = const.BAR_SIO2;
					obj.mass = const.MASS_SIO2;
					obj.massxy = const.MASS_SIO2;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 0;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 17;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'nc-Si')
					obj.name = "nc-Si.";
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = const.MASS_iSI;
					obj.massxy = const.MASS_AL;
					obj.valley = 1;
					obj.cond = 1;
					obj.Q = Q;
					obj.Ef = -const.EG_SI/2+const.KB*const.TEMP*log((obj.Q*1e6)/const.NI)/const.ELEC;
					obj.base = const.BASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 18;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'nSi100Z')
					obj.name = "nSi100Z";
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = const.MASS_SI_Z;
					obj.massxy = const.MASS_SI_XY;
					obj.valley = const.NUMVALLY_SI_Z;
					obj.cond = 1;
					obj.Q = Q;
					obj.Ef = -const.EG_SI/2+const.KB*const.TEMP*log((obj.Q*1e6)/const.NI)/const.ELEC;
					obj.base = const.BASE;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 19;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'nSi100XY')
					obj.name = "nSi100XY";
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = const.MASS_SI_XY;
					obj.massxy = const.MASS_SI_Z;
					obj.valley = const.NUMVALLY_SI_XY;
					obj.cond = 1;
					obj.Q = Q;
					obj.Ef = -const.EG_SI/2+const.KB*const.TEMP*log((obj.Q*1e6)/const.NI)/const.ELEC;
					obj.base = BASE;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 20;
                    obj.d = ML*const.ML;
				% elseif strcmp(materialName, 'n-Si.Sub.')
				% 	obj.name="n-Si.Sub.";
				% 	obj.die=DIE_iSI  * DIEELECSTAR;
				% 	obj.bar=BAR_iSI;
				% 	obj.mass=MASS_iSI;
				% 	obj.massxy=MASS_AL;
				% 	obj.valley=NUMVALLY_SI;
				% 	obj.cond=0;
				%	obj.smt = 21;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'CaF2(4ML)')
					obj.name = "CaF2(4ML)";
					obj.die = const.DIE_CAF2 * const.DIEELECSTAR;
					obj.bar = 1.6;
					obj.mass = 0.85;
					obj.massxy = 0.85;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 0;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 22;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, 'p-i-Si')
					obj.name = "p-i-Si";
					obj.die = const.DIE_iSI  * const.DIEELECSTAR;
					obj.bar = const.BAR_iSI;
					obj.mass = 0.55;
					obj.massxy = const.MASS_iSI;
					obj.valley = const.NUMVALLY_SI;
					obj.cond = 0;	
					obj.Q = Q;
					nv = 2*pow(const.MSTAR*const.KB*const.TEMP/const.HBAR/const.HBAR/pi/2.0, 1.5 ) * (pow(0.16,1.5) + pow(0.49,1.5));
					obj.Ef = Ed + const.KB*const.TEMP*log(-1.0/4.0 + pow(1+8*obj.Q*1e6*exp(-Ed/const.KB/const.TEMP)/nv, 0.5)/4)/const.ELEC;
					obj.base = const.pBASE;
					obj.Work = -obj.Ef - obj.bar + obj.base;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 23;
                    obj.d = ML*const.ML;
				elseif strcmp(materialName, '3MLCaF2')
					obj.name = "3MLCaF2";
					obj.die = const.DIE_CAF2 * const.DIEELECSTAR;
					obj.bar = 1.5;
					obj.mass = 1.0;
					obj.massxy = 1.0;
					obj.valley = 1;
					obj.cond = 3;
					obj.Ef = 0;
					obj.Q = Q;
					obj.base = const.pBASE;
					obj.Work = 0;
					obj.divnum = ML*const.DX;
					obj.NX = NX;
					obj.smt = 24;
                    obj.d = ML*const.ML;
				end
			end
		end
	end
end
