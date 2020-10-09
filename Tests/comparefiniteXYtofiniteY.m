% close all
closewaitbars
clearvars -except array*

f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;
z0 = Constants.z0;

f = 24 * 1e9;

th = eps * pi/180;
ph = 0 * pi/180;

dx = 0.45*l0;
dy = dx;
% wslot = 0.05*l0;
% dslot = 0.05*l0;
wslot = 1.1e-3;
dslot = 2e-3;
walled = 0;
dedge = 0.25*l0;
zfeed = 100;


Nx = 10;
Ny = 10;
excitation = ones(Nx, Ny);

tlineup = FreeSpace();
tlinedown = FreeSpace();

% Backing reflector
erback = 2.2;
hback = 1.9e-3;
tlinedown = ShortedLine(erback, hback, erback);

% ESA ADL
z1 = 90;
z2 = z0;
p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 19e9;
f0design = 29e9;
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);
tlineup = TerminatedTLine(slab, FreeSpace());


slot = Slot(dx, dy, wslot, dslot, walled);

arrayInf = InfiniteArray(slot, tlineup, tlinedown);

%% Compare D
% FiniteArray
if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
    array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
end
array.InitializeDs(f);

% FiniteArrayY
if(~exist('arrayY', 'var') || ~isa(arrayY, 'FiniteArrayY') || arrayY.Ny ~= Ny)
    arrayY = FiniteArrayY(slot, tlineup, tlinedown, Ny, excitation(1,:), zfeed);
end
arrayY.InitializeDs(f, th, ph);

% if(~exist('array', 'var'))
%     array = load('H:\Git\PhD-Matlab\Validations\finitearray_ESA\10x10.mat');
%     array = array.array;
%     arrayY = load('H:\Git\PhD-Matlab\Temp\arrayY_10ESA.mat');
%     arrayY = arrayY.arrayY;
% end

mxs = [-arrayY.numM:arrayY.numM];
nyp = 0;

fiY = arrayY.Dmy_fs == f;
fi = array.D_fs == f;

th = arrayY.Dmy_th;
ph = arrayY.Dmy_ph;

dx = arrayY.unitcell.dx;

[k0, kx0, ~, ~] = k(f, 1, th, ph);
kxm = kx0 - 2*pi*mxs/dx;

kxs = linspace(min(kxm), max(kxm), 1000);
kxs = sort([kxs kxm]);

Dmys = arrayY.Dmys(:, 1, fiY);

deformedpath = 0; % Nondeformed path.
Ds = array.D_interpolants{deformedpath+1, fi, nyp+1}(kxs).';

[hFig, hAx] = figureex;
    repeatcolormap(hAx, 2);
    plot(real(kxm)./k0, real(Dmys), 'x');
    plot(real(kxm)./k0, imag(Dmys), 'o');
    
    plot(real(kxs)./k0, real(Ds));
    plot(real(kxs)./k0, imag(Ds), '--');
    title('D');
    
Dsm = array.D_interpolants{deformedpath+1, fi, nyp+1}(kxm).';
[hFig, hAx] = figureex;
    repeatcolormap(hAx, 2);
    plot(real(kxm)./k0, 100.*real(Dmys - Dsm)./real(Dsm), '');
    plot(real(kxm)./k0, 100.*imag(Dmys - Dsm)./imag(Dsm), '--');
    title('D: Difference in %');
    
%% Compare Ky Integral (Gfa)
arrayY.InitializeZMatrix(f, th, ph);
array.InitializeKyInts(f);

fiY = arrayY.ints_fs == f;
fi = array.KyInt_fs == f;

deformedpath = 0; % Nondeformed path.
dny = 0;
kxs = array.KyInt_kxs{deformedpath+1, fi, dny+1};
kyints_fY = arrayY.ints(:, 1, fiY);
kyints_f = array.KyInt_interpolants{deformedpath+1, fi, dny+1}(kxs).';

[hFig, hAx] = figureex;
    repeatcolormap(hAx, 2);
    plot(real(kxm)./k0, real(kyints_fY), 'x');
    plot(real(kxm)./k0, imag(kyints_fY), 'o');
    
    plot(real(kxs)./k0, real(kyints_f));
    plot(real(kxs)./k0, imag(kyints_f), '--');
    title('G^{FA}');
    
kyints_fm = array.KyInt_interpolants{deformedpath+1, fi, dny+1}(kxm).';
[hFig, hAx] = figureex;
    repeatcolormap(hAx, 2);
    plot(real(kxm)./k0, 100.*real(kyints_fY - kyints_fm)./real(kyints_fm), '');
    plot(real(kxm)./k0, 100.*imag(kyints_fY - kyints_fm)./imag(kyints_fm), '--');
    title('G^{FA}: Difference in %');

%% Compare Z matrix
fs = (20:2:30) * 1e9;
fs = [fs, 12e9, 14e9, 16e9, 18e9, 32e9];
fs = sort(fs);

arrayY.InitializeZMatrix(fs, th, ph);
array.InitializeZMatrix(fs);

ZmatY = arrayY.Zmat;
Zmat = array.ReducedZMatrix(fs);

Z1xY = ZmatY(1,:,:);
nx = floor(Nx/2);
ny = 0;
nyp = 0:Ny-1;
nxp = floor(Nx/2);
for(fi = 1:length(fs))
    Z1x{fi} = Zmat{fi}((nx+1)+Nx*ny, (nxp+1)+Nx*nyp);
end
Z1x = reshape(cell2mat(Z1x), Ny, length(fs));


%% Compare active Z
ZasInf = arrayInf.GetInputImpedance(fs, th, ph);
Zas = array.GetInputImpedance(fs, excitation);
ZasY = arrayY.GetInputImpedance(fs, th, ph);

[hFig, hAx] = figureex;
    repeatcolormap(hAx, 2);
    for(nx = floor(Nx/2):floor(Nx/2)+1)
        for(ny = 1:Ny)
            hPlot = plot(hAx, fs./1e9, real(squeeze(Zas(nx,ny,:))));
            addlegendentry(hAx, hPlot, sprintf('%i, %i', nx, ny));
            plot(hAx, fs./1e9, imag(squeeze(Zas(nx,ny,:))), '--');
        end
    end
    for(ny = 1:Ny)
        hPlot = plot(hAx, fs./1e9, real(squeeze(ZasY(ny,:))), 'x');
        addlegendentry(hAx, hPlot, sprintf('%i', ny));
        plot(hAx, fs./1e9, imag(squeeze(ZasY(ny,:))), 'o');
    end
    hPlot = plot(hAx, fs./1e9, ZasInf, 'k');
    addlegendentry(hAx, hPlot, 'Infinite');
    plot(hAx, fs./1e9, imag(squeeze(ZasInf)), 'k:');
    
% [hFig, hAx] = figureex;
%     repeatcolormap(hAx, 2);
%     plot(fs, 100.*real(Zas - ZasY)./real(Zas), '');
%     plot(fs, 100.*imag(Zas - ZasY)./imag(Zas), '--');
%     title('Difference in %');
    
%% Compare active S
SInf = (ZasInf - zfeed) ./ (ZasInf + zfeed);
S = (Zas - zfeed) ./ (Zas + zfeed);
SY = (ZasY - zfeed) ./ (ZasY + zfeed);

[hFig, hAx] = figureex;
    for(nx = floor(Nx/2):floor(Nx/2)+1)
        for(ny = 1:Ny)
            hPlot = plot(hAx, fs./1e9, 20*log10(abs(squeeze(S(nx,ny,:)))));
            addlegendentry(hAx, hPlot, sprintf('%i, %i', nx, ny));
        end
    end
    for(ny = 1:Ny)
        hPlot = plot(hAx, fs./1e9, 20*log10(abs(squeeze(SY(ny,:)))), 'x');
        addlegendentry(hAx, hPlot, sprintf('%i', ny));
    end
    hPlot = plot(hAx, fs./1e9, 20*log10(abs(SInf)), 'k');
    addlegendentry(hAx, hPlot, 'Infinite');

    

%%
% Dmy = [0.0200446886612149 - 362.073166635821i,0.0106033914525551 + 0.714074027550655i;0.0183120409174059 - 346.125953497723i,0.00968684251769413 + 0.652404484101316i;0.0166307370424242 - 330.104924386618i,0.00879745381239659 + 0.592553777110509i;0.0150048953273331 - 314.003323262028i,0.00793740386290671 + 0.534668990705014i;0.0134387773236505 - 297.812803916219i,0.00710894697894285 + 0.478902424354931i;0.0119367641568105 - 281.522622508649i,0.00631440072375877 + 0.425410750474708i;0.0105033277677378 - 265.118517401839i,0.00555613070132549 + 0.374353988354760i;0.00914299699160985 - 248.581211348249i,0.00483653261239388 + 0.325894290740466i;0.00786031854375864 - 231.884482900108i,0.00415801161643751 + 0.280194545189052i;0.00665981317199315 - 214.992776550220i,0.00352295913681218 + 0.237416799206218i;0.00554592744731098 - 197.858343406188i,0.00293372735848680 + 0.197720525906118i;0.00452298189284274 - 180.417889374948i,0.00239260178891323 + 0.161260755259234i;0.00359511638291586 - 162.588574710994i,0.00190177237462169 + 0.128186104489437i;0.00276623396662345 - 144.262827663162i,0.00146330378454290 + 0.0986367493533868i;0.00203994446765449 - 125.300710039513i,0.00107910557484993 + 0.0727423853087830i;0.00141950936804481 - 105.517841809535i,0.000750903033078718 + 0.0506202333485776i;0.000907789582161133 - 84.6684960390805i,0.000480209551075775 + 0.0323731481255900i;0.000507197756766698 - 62.4388992148854i,0.000268301391384890 + 0.0180876360501551i;0.000219656668697919 - 38.5570140422507i,0.000116195684041850 + 0.00776283906727567i;4.65651971525989e-05 - 13.6286980633601i,2.46324183062784e-05 - 0.0144436536209983i;2.91629895909490 + 4.06292565556162i,-0.667618077157650 - 0.774146846229539i;4.65651971525989e-05 - 13.6286980633601i,2.46324183062784e-05 - 0.0144436536209983i;0.000219656668697919 - 38.5570140422507i,0.000116195684041850 + 0.00776283906727567i;0.000507197756766698 - 62.4388992148854i,0.000268301391384890 + 0.0180876360501551i;0.000907789582161133 - 84.6684960390805i,0.000480209551075775 + 0.0323731481255900i;0.00141950936804481 - 105.517841809535i,0.000750903033078718 + 0.0506202333485776i;0.00203994446765449 - 125.300710039513i,0.00107910557484993 + 0.0727423853087830i;0.00276623396662345 - 144.262827663162i,0.00146330378454290 + 0.0986367493533868i;0.00359511638291586 - 162.588574710994i,0.00190177237462169 + 0.128186104489437i;0.00452298189284274 - 180.417889374948i,0.00239260178891323 + 0.161260755259234i;0.00554592744731098 - 197.858343406188i,0.00293372735848680 + 0.197720525906118i;0.00665981317199315 - 214.992776550220i,0.00352295913681218 + 0.237416799206218i;0.00786031854375864 - 231.884482900108i,0.00415801161643751 + 0.280194545189052i;0.00914299699160985 - 248.581211348249i,0.00483653261239388 + 0.325894290740466i;0.0105033277677378 - 265.118517401839i,0.00555613070132549 + 0.374353988354760i;0.0119367641568105 - 281.522622508649i,0.00631440072375877 + 0.425410750474708i;0.0134387773236505 - 297.812803916219i,0.00710894697894285 + 0.478902424354931i;0.0150048953273331 - 314.003323262028i,0.00793740386290671 + 0.534668990705014i;0.0166307370424242 - 330.104924386618i,0.00879745381239659 + 0.592553777110509i;0.0183120409174059 - 346.125953497723i,0.00968684251769413 + 0.652404484101316i;0.0200446886612149 - 362.073166635821i,0.0106033914525551 + 0.714074027550655i];
% kxm = [28856.1102996396;27413.3047846576;25970.4992696756;24527.6937546936;23084.8882397117;21642.0827247297;20199.2772097477;18756.4716947657;17313.6661797837;15870.8606648018;14428.0551498198;12985.2496348378;11542.4441198558;10099.6386048739;8656.83308989187;7214.02757490990;5771.22205992792;4328.41654494594;2885.61102996396;1442.80551498198;1.94799325333888e-15;-1442.80551498198;-2885.61102996396;-4328.41654494594;-5771.22205992792;-7214.02757490990;-8656.83308989187;-10099.6386048739;-11542.4441198558;-12985.2496348378;-14428.0551498198;-15870.8606648018;-17313.6661797837;-18756.4716947657;-20199.2772097477;-21642.0827247297;-23084.8882397117;-24527.6937546936;-25970.4992696756;-27413.3047846576;-28856.1102996396];
% D = [0.0754815959567781 - 362.072560537287i,0.0103857654796190 + 0.714092059148305i;0.0739951479670270 - 346.125884234002i,0.00947510205309267 + 0.652409191814338i;0.0725825971789679 - 330.104305616897i,0.00859250697858434 + 0.592571821590970i;0.0712514512523598 - 314.002163212869i,0.00773975938355220 + 0.534699529837460i;0.0700125895576317 - 297.811271978767i,0.00691906368685588 + 0.478940209547905i;0.0688804031301978 - 281.520954536541i,0.00613272877440260 + 0.425449041976185i;0.0678738134832565 - 265.117076711181i,0.00538310669675681 + 0.374384482240454i;0.0670176080375047 - 248.580602899209i,0.00467257472811586 + 0.325907088725010i;0.0663465340407973 - 231.883481874359i,0.00400402559418243 + 0.280212070596063i;0.0659025224838621 - 214.989996186709i,0.00337974186874509 + 0.237454727650724i;0.0657363292754165 - 197.853806382213i,0.00280170730466122 + 0.197770888928210i;0.0659124901913163 - 180.411987579374i,0.00227216542020397 + 0.161314009767463i;0.0665061302005270 - 162.586684725245i,0.00179325897457040 + 0.128231243509000i;0.0676286813097795 - 144.261067535690i,0.00136700116321120 + 0.0986614832083697i;0.0694094541972877 - 125.299757524839i,0.000995538886035694 + 0.0727529317026913i;0.0720353398858096 - 105.517185489418i,0.000681098685355433 + 0.0506595615444606i;0.0757521091400071 - 84.6671955538835i,0.000424334482784421 + 0.0324303115081307i;0.0807436075785565 - 62.4386341344158i,0.000226440308889436 + 0.0181504544464178i;0.0862152944690705 - 38.5559210655159i,8.62532635389839e-05 + 0.00777881417672671i;0.0853025208338620 - 13.6293063070759i,-0.000313746571038852 - 0.0161140674172170i;2.91665468497095 + 4.06136161497206i,-0.667262478392299 - 0.773947424544293i;0.0853025208338620 - 13.6293063070759i,-0.000313746571038664 - 0.0161140674172165i;0.0862152944690706 - 38.5559210655159i,8.62532635395473e-05 + 0.00777881417672783i;0.0807436075785565 - 62.4386341344158i,0.000226440308888323 + 0.0181504544464185i;0.0757521091400070 - 84.6671955538836i,0.000424334482786975 + 0.0324303115081329i;0.0720353398858096 - 105.517185489418i,0.000681098685354865 + 0.0506595615444622i;0.0694094541972877 - 125.299757524839i,0.000995538886035796 + 0.0727529317026905i;0.0676286813097795 - 144.261067535690i,0.00136700116321296 + 0.0986614832083677i;0.0665061302005271 - 162.586684725245i,0.00179325897456577 + 0.128231243509002i;0.0659124901913161 - 180.411987579374i,0.00227216542020284 + 0.161314009767461i;0.0657363292754167 - 197.853806382213i,0.00280170730466147 + 0.197770888928214i;0.0659025224838623 - 214.989996186710i,0.00337974186874556 + 0.237454727650728i;0.0663465340407973 - 231.883481874359i,0.00400402559418242 + 0.280212070596062i;0.0670176080375049 - 248.580602899209i,0.00467257472811722 + 0.325907088725004i;0.0678738134832564 - 265.117076711181i,0.00538310669675118 + 0.374384482240450i;0.0688804031301976 - 281.520954536541i,0.00613272877439530 + 0.425449041976180i;0.0700125895576320 - 297.811271978767i,0.00691906368685551 + 0.478940209547907i;0.0712514512523596 - 314.002163212869i,0.00773975938355632 + 0.534699529837463i;0.0725825971789680 - 330.104305616897i,0.00859250697857148 + 0.592571821590968i;0.0739951479670266 - 346.125884234002i,0.00947510205309209 + 0.652409191814341i;0.0754815959567788 - 362.072560537287i,0.0103857654796142 + 0.714092059148284i];
% kxs = [28856.1102996396;27413.3047846576;25970.4992696756;24527.6937546936;23084.8882397117;21642.0827247297;20199.2772097477;18756.4716947657;17313.6661797837;15870.8606648018;14428.0551498198;12985.2496348378;11542.4441198558;10099.6386048739;8656.83308989187;7214.02757490990;5771.22205992792;4328.41654494594;2885.61102996396;1442.80551498198;1.94799325333888e-15;-1442.80551498198;-2885.61102996396;-4328.41654494594;-5771.22205992792;-7214.02757490990;-8656.83308989187;-10099.6386048739;-11542.4441198558;-12985.2496348378;-14428.0551498198;-15870.8606648018;-17313.6661797837;-18756.4716947657;-20199.2772097477;-21642.0827247297;-23084.8882397117;-24527.6937546936;-25970.4992696756;-27413.3047846576;-28856.1102996396];
% 
% [hFig, hAx] = figureex;
%     repeatcolormap(hAx, 2);
%     plot(real(kxm)./k0, real(Dmy(:,1)), 'x');
%     plot(real(kxm)./k0, imag(Dmy(:,1)), 'o');
%     
%     plot(real(kxs)./k0, real(D(:,1)));
%     plot(real(kxs)./k0, imag(D(:,1)), '--');
