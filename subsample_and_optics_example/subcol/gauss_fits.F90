module mod_gauss_fits
   implicit none
contains

!add time info as an input eventually. Also, add do_ls_param, ls_p_rate
subroutine gauss_fits(npoints,nlev,ncol,cv_p_rate,do_cv_param,do_max_p_in_col, &
                      millisecond,second,minute,cv_col_param) 

USE COSP_KINDS, ONLY: wp
      INTEGER :: npoints,    &    ! Number of model points in the horizontal
                 nlev,       &    ! Number of model levels in column
                 ncol             ! Number of subcolumns per model grid point
      REAL(WP),dimension(100) :: cdf, frac_bins, z_eqn, yfit_eqn, yfit_sum
      REAL(WP),dimension(120) :: rain_bins, coeff1, coeff2, coeff3, coeff4,coeff5,coeff6 !Array
      INTEGER, dimension(1)   :: coeff_bin, rand_num_ind
      INTEGER                 :: i,j,t
      REAL(WP)                :: millisecond, second, minute!, cloud_frac_val !scalar value !**Eventually make these inputs
      INTEGER,dimension(:),allocatable :: seedB
      INTEGER,dimension(npoints)       :: cv_nearsfc_loc !array location of near surface rain
                                          !ls_nearsfc_loc
      INTEGER,dimension(npoints)       :: cv_col_param !,ls_frac !number of subcols that should be filled 
      REAL(WP),dimension(npoints)      :: rand_num_for_cldfrac, & !random numbers to select cloud fraction based on cdf
                                          cv_frac,              & ! cloud fraction
                                          col_cv_prec!,        !& !Either max or near sfc rain rate in each column
                                          !ls_frac,             &
                                          !col_ls_prec
      REAL(WP),dimension(npoints,nlev) :: cv_p_rate, & 
                                          cv_mm_day, &  !used for cv/ls_col_frac param
                                          temp_cv_mm_day !Used for reversing level to go from sfc to TOA
                                          !ls_p_rate
                                          !ls_mm_day    !used for cv/ls_col_frac param
      LOGICAL :: do_cv_param, do_max_p_in_col

print *, 'millisecond', millisecond
print *, 'second', second
print *, 'minute', minute
!Fill arrays of cloud fracion and rain rates that are used to determine each CDF per rain rate bin
      frac_bins = (/ (I, I = 1, 100) /)
      frac_bins = (frac_bins/100.) - .005

      rain_bins = (/ (J, J = 1, 120)/)
      rain_bins = rain_bins - 0.5

!!Determine precipitaiton rate to use at each point. Use either the max rain rate in each column or the near surface rain
      if (do_cv_param) then

        cv_mm_day = cv_p_rate*86400. !put rainrate in mm/day, since TRMM fit based on mm/day

        if (do_max_p_in_col) then
         col_cv_prec = maxval(cv_mm_day, dim = 2) !fine max precip rate per column
        else !Default is to find precip closest to surface n each column 

         !reverse the levels of the rain rate so they go surfac to TOA. That way when we mask the array to find the
         !rain closest to the surface the maxloc will give the correct array index b/c maxloc gives the first instance
         !of where the max value is. Currently the cv/ls_mm_day variables are TOA to surface. ERileyDellaripa (2/4/2019) 
         temp_cv_mm_day = cv_mm_day(:,nlev:1:-1)
         cv_mm_day      = cv_mm_day(:,nlev:1:-1)
         where(temp_cv_mm_day .gt. 0) temp_cv_mm_day = 1
         cv_nearsfc_loc = maxloc(temp_cv_mm_day, dim = 2)

         do j = 1,npoints
          col_cv_prec(j) = cv_mm_day(j, cv_nearsfc_loc(j))
         enddo
        endif !if(do_max_p_in_col) 
print *, 'col_cv_prec(1:20)=', col_cv_prec(1:20)
!Gaussian Coefficients
!have an if statement to select convective coefficients vs. large-scale/stratiform coefficients
!Below is just coefficients for convective over land and ocean combined
        coeff1 = (/705.873,      23.5916,      15.8820,      12.4454,      10.6658,      9.77671, &
      9.14880,      8.77748,      8.68286,      8.50575,      8.40796,      8.45044, &
      8.42248,      8.29453,      8.28179,      8.53937,      8.45761,      8.27237, &
      8.45939,      8.38701,      8.32310,      8.46185,      8.43146,      8.28767, &
      8.56599,      8.36547,      8.44702,      8.42221,      8.08302,      8.15192, &
      8.31430,      8.22149,      7.94934,      7.89766,      7.98632,      8.07544, &
      8.25443,      7.79529,      7.86957,      7.97676,      8.01271,      8.17470, &
      7.78040,      7.86477,      7.65418,      7.83410,      7.75081,      7.79590, &
      7.53953,      7.57709,      7.82730,      7.77789,      7.60649,      7.87077, &
      7.69078,      7.11822,      7.92643,      7.82562,      7.85921,      7.34736, &
      6.97912,      7.39185,      7.41789,      7.59665,      7.83264,      7.09878, &
      7.27298,      7.26285,      7.40040,      7.06468,      7.24910,      6.91211, &
      6.86139,      6.94590,      6.93427,      7.36334,      8.51558,      7.53109, &
      6.94756,      7.80104,      7.11426,      7.09944,      7.28409,      7.17016, &
      7.14061,      7.04041,      6.39235,      6.75556,      6.73471,      6.56398, &
      6.66908,      7.41178,      7.01170,      7.52170,      6.61242,      6.88643, &
      6.76742,      7.33465,      6.95158,      8.70462,      6.71506,      8.19908, &
      9.45001,      9.29065,      7.19872,      8.35447,      7.44559,      5.89912, &
      7.78527,      6.22007,      6.83044,      6.09346,      13.6463,      20.6898, &
      5.89651,      6.62073,      6.48200,      6.52083,      18.1614,      12.0642/)
              
        coeff2=(/-0.0382389,    0.0414715,    0.0602159,    0.0702888,    0.0767319,    0.0801908, &
    0.0817204,    0.0843448,    0.0853339,    0.0861372,    0.0875907,    0.0908833, &
    0.0896298,    0.0926564,    0.0934401,    0.0941976,    0.0968262,    0.0941161, &
    0.0959174,    0.0982337,     0.100369,     0.101899,     0.103374,     0.105277, &
     0.105526,     0.103072,     0.105157,     0.108972,     0.107593,     0.109813, &
     0.110596,     0.112515,     0.112903,     0.114051,     0.114871,     0.116331, &
     0.118577,     0.117299,     0.115912,     0.121174,     0.119837,     0.121859, &
     0.119031,     0.118996,     0.124054,     0.126346,     0.123883,     0.127538, &
     0.128258,     0.120017,     0.125854,     0.135763,     0.131017,     0.135069, &
     0.134774,     0.136277,     0.125999,     0.133921,     0.130122,     0.137161, &
     0.141159,     0.135028,     0.143356,     0.140844,     0.138165,     0.135608, &
     0.146621,     0.140035,     0.147830,     0.143024,     0.137535,     0.148335, &
     0.150574,     0.148852,     0.146740,     0.148268,     0.139818,     0.149635, &
     0.147821,     0.150314,     0.155333,     0.148683,     0.139687,     0.157440, &
     0.146726,     0.149897,     0.166739,     0.163705,     0.159137,     0.161381, &
     0.160399,     0.171204,     0.151055,     0.151714,     0.127625,     0.163433, &
     0.167350,     0.162127,     0.176422,     0.174594,     0.173338,     0.116096, &
     0.203015,     0.153569,     0.168835,     0.145302,     0.167564,     0.182136, &
     0.159155,     0.178209,     0.182447,     0.178232,     0.136400,     0.111742, &
     0.183406,     0.173126,     0.171569,     0.179352,     0.126200,     0.237941/)

        coeff3=(/0.0199818,    0.0169630,    0.0263397,    0.0351587,    0.0417259,    0.0447758, &
    0.0449863,    0.0459215,    0.0449362,    0.0447091,    0.0457754,    0.0453509, &
    0.0456923,    0.0462062,    0.0463306,    0.0441507,    0.0464822,    0.0460183, &
    0.0457673,    0.0470492,    0.0484009,    0.0478511,    0.0466429,    0.0493174, &
    0.0457618,    0.0464918,    0.0481797,    0.0498743,    0.0511409,    0.0517723, &
    0.0493492,    0.0469103,    0.0515994,    0.0544480,    0.0511178,    0.0526332, &
    0.0485076,    0.0532411,    0.0523258,    0.0490911,    0.0521942,    0.0526147, &
    0.0523247,    0.0549065,    0.0635183,    0.0548231,    0.0554595,    0.0574905, &
    0.0594238,    0.0527295,    0.0567083,    0.0584718,    0.0527863,    0.0516353, &
    0.0555525,    0.0548879,    0.0520969,    0.0527599,    0.0507838,    0.0634184, &
    0.0656692,    0.0598679,    0.0595362,    0.0545073,    0.0650621,    0.0594522, &
    0.0602917,    0.0603958,    0.0644794,    0.0667981,    0.0568071,    0.0740435, &
    0.0651259,    0.0677056,    0.0612082,    0.0683696,    0.0494006,    0.0660586, &
    0.0706698,    0.0485437,    0.0621373,    0.0655728,    0.0568298,    0.0609397, &
    0.0606096,    0.0564775,    0.0708071,    0.0601342,    0.0685513,    0.0720213, &
    0.0706667,    0.0483823,    0.0600589,    0.0578166,   0.00974120,    0.0608004, &
    0.0734310,    0.0600506,   0.00942726,    0.0208097,    0.0512615,   0.00903693, &
    0.0135501,    0.0340707,    0.0590975,    0.0527556,    0.0593223,    0.0561646, &
    0.0463545,    0.0687764,    0.0589320,    0.0619222,    0.0122567,   0.00370934, &
    0.0632564,    0.0607365,    0.0716710,    0.0708045,   0.00127102,   0.00365503/)


        coeff4 =(/0.467342,   0.00529722,    -0.240571,    -0.482603,    -0.550525,    -0.462593, &
    -0.161937,   -0.0714611,    0.0541969,     0.152619,    0.0862023,    0.0929808, &
    0.0407834,    0.0715796,    0.0348292,    0.113093,   -0.0876833,    0.0607330, &
   -0.0282712,    -0.146555,    -0.216726,    -0.267125,    -0.116516,    -0.290908, &
   -0.0950166,    -0.104707,    -0.317058,    -0.435161,    -0.397234,    -0.480632, &
    -0.354192,   -0.0570672,    -0.348671,    -0.589398,    -0.357915,    -0.528159, &
    -0.259965,    -0.409887,    -0.411962,   -0.0920201,    -0.501202,    -0.591189, &
    -0.332160,    -0.663124,     -1.19168,    -0.526741,    -0.557786,    -0.780238, &
    -0.733800,    -0.314603,    -0.802176,    -0.832842,    -0.323522,    -0.354523, &
    -0.538084,    -0.156537,    -0.512650,    -0.386808,    -0.355908,     -1.04006, &
    -0.981300,    -0.759345,    -0.699840,    -0.487163,     -1.49676,    -0.632387, &
    -0.660282,    -0.758948,     -1.10043,     -1.07690,    -0.475742,     -1.50597, &
    -0.828725,    -0.972613,    -0.567546,     -1.35483,    -0.464202,     -1.34214, &
     -1.29205,   -0.0954215,    -0.730005,     -1.03518,    -0.551488,    -0.799317, &
    -0.794264,    -0.452372,    -0.900877,    -0.230781,    -0.871624,     -1.02715, &
    -0.933969,     0.194633,    -0.576410,    -0.698817,      3.04112,    -0.415142, &
     -1.18034,    -0.680452,      3.09745,      1.87184,     0.157241,     2.56647, &
      2.39946,     0.514038,    -0.448912,    -0.702939,    -0.674713,    0.0662954, &
    -0.168376,    -0.626720,     0.117288,    -0.232882,      1.86109,      2.56609, &
    -0.158036,    -0.297233,    -0.899888,    -0.922341,      2.42741,      2.94993/)

        coeff5 =(/-1.78015,    0.0117685,     0.840815,      1.71920,      1.99710,    1.77656, &
     0.878911,     0.613223,     0.252118,   -0.0201119,     0.211463,     0.160972, &
     0.410082,     0.264405,     0.439477,     0.154211,     0.787061,     0.365947, &
     0.629388,      1.04953,     1.20289,      1.38184,     0.876462,      1.40190, &
     0.771855,     0.929478,      1.56413,      1.82923,      1.83719,      2.03965, &
      1.63951,     0.714865,      1.63733,      2.39975,      1.70905,      2.16193, &
      1.36197,      1.80947,     1.91139,     0.727172,      2.16853,      2.32092, &
      1.59752,      2.69766,     4.13571,      1.99980,      2.17374,      2.85295, &
      2.64416,      1.69318,     3.07515,      2.91573,      1.63816,      1.64742, &
      2.08727,      1.06860,     2.32824,      1.66713,      1.85504,      3.74307, &
      3.63461,      2.86040,     2.53545,      2.13049,      4.98668,      2.72505, &
      2.41882,      2.94014,     3.71553,      3.78402,      2.16010,      5.04479, &
      3.07629,      3.27195,     2.31112,      4.39917,      1.82003,      4.43017, &
      4.41013,     0.998325,      2.66833,      3.67222,      2.47454,      3.11376, &
      3.23254,      2.34634,      3.28331,     0.995452,      2.92802,      3.46116, &
      3.01347,    0.0254561,      2.39485,      2.66926,     -6.74949,      1.62806, &
      3.63369,      2.38462,     -7.10985,     -3.88679,     0.544872,     -4.60702, &
     -5.18879,    -0.525593,      1.59615,      2.52306,      2.34070,      1.27076, &
      1.76677,      2.45559,    -0.665180,      1.69917,     -3.47350,     -4.61712, &
      1.41691,      1.40866,      2.97881,      3.11495,     -3.20854,     -5.82279/)

        coeff6 = (/1.44669,   -0.0231499,    -0.651579,     -1.34764,     -1.58156,     -1.44934, &
    -0.821222,    -0.636956,    -0.391187,    -0.210907,    -0.390173,    -0.338006,  &
    -0.560054,    -0.429273,    -0.586624,    -0.359667,    -0.819333,    -0.535930,  &
    -0.718704,     -1.04793,     -1.12635,     -1.26621,    -0.883487,     -1.25381,  &
    -0.789751,    -0.966295,     -1.41142,     -1.55193,     -1.61808,     -1.73797,  &
     -1.44395,    -0.777748,     -1.45071,     -2.00702,     -1.52236,     -1.81173,  &
     -1.25250,     -1.56285,     -1.68675,    -0.739662,     -1.86040,     -1.90930,  &
     -1.42648,     -2.25578,     -3.19219,     -1.61530,     -1.77888,     -2.26472,  &
     -2.08229,     -1.56644,     -2.50030,     -2.26072,     -1.48656,     -1.45031,  &
     -1.70520,     -1.05239,     -2.03577,     -1.42801,     -1.70080,     -2.94519,  &
     -2.90038,     -2.30065,     -1.99846,     -1.83303,     -3.75776,     -2.32863,  &
     -1.92033,     -2.39846,     -2.82136,     -2.93712,     -1.88597,     -3.81137,  &
     -2.45392,     -2.47748,     -1.92964,     -3.26333,     -1.49366,     -3.31628,  &
     -3.36424,     -1.05623,     -2.11626,     -2.86293,     -2.14876,     -2.54261,  &
     -2.69561,     -2.14079,     -2.59466,    -0.845992,     -2.21376,     -2.62103,  &
     -2.22139,    -0.318227,     -2.01111,     -2.16367,      3.51654,     -1.33404,  &
     -2.60435,     -1.83641,      3.87938,      1.85239,    -0.876665,      1.65398,  &
      2.62142,    -0.134083,     -1.23991,     -1.98053,     -1.80063,     -1.59510,  &
     -1.85812,     -2.00498,     0.621876,     -1.68194,      1.36920,      1.67800,  &
     -1.45093,     -1.23868,     -2.23285,     -2.35794,     0.182382,      2.55784/)

        !1) Generate Random number to select point along CDF
        call random_seed(size=t)
        !print *, 't = ', t

        ALLOCATE(seedB(1:t))

        seedB(t)=millisecond
        seedB(1)=millisecond*second*minute

        !print *, 'seedB=', seedB

        call random_seed(put = seedB)

        deallocate(seedB)

        call random_number(rand_num_for_cldfrac)
        print *, 'random_number(1:15) = ', rand_num_for_cldfrac(1:15)

        !**************************************************************************************************
        !Time to loop through Npoints
        !**************************************************************************************************
        do j=1,npoints
          !2) Determine which rain_bin the CDF will need to be generated for by locating the bin the 
          !rainrate falls in. That bin will be used to select the bin (i.e., array element) in the Gaussian 
          !coefficient array
          coeff_bin = minloc(abs(col_cv_prec(j) - rain_bins))

          !3) Determine PDF using Gaussian coefficients
          !need to use minval to make coeffX(coeff_bin) a scalar instead of an array of rank & shape 1
          z_eqn    = (frac_bins - minval(coeff2(coeff_bin)))/minval(coeff3(coeff_bin))
          yfit_eqn = minval(coeff1(coeff_bin))*exp((0-z_eqn**2)/2.) + minval(coeff4(coeff_bin)) + &
                     (minval(coeff5(coeff_bin))*frac_bins) + (minval(coeff6(coeff_bin))*frac_bins**2)

          !4) Compute CDF
          DO i=1,size(yfit_eqn)
            yfit_sum(i) = sum(yfit_eqn(1:i))
          ENDDO

          !print *, 'yfit_sum', yfit_sum

          CDF = yfit_sum/sum(yfit_eqn)

          !5) Isolate location of random number in CDF
           rand_num_ind = minloc(abs(rand_num_for_cldfrac(j) - cdf))

          !6) Select cloud fraction associated with the fraction value from the random number
          cv_frac(j) = minval(frac_bins(rand_num_ind))

        enddo !loop over Npoints to find cloud fraction

        !7) Determine number of subcolumns that should be filled based on cv_frac
        cv_col_param = NINT(cv_frac*ncol)
     
        !ensures 1 col has prec if cv prec 
        where(cv_col_param .lt. 1 .and. col_cv_prec .gt. 0) cv_col_param = 1 
        !ensure cv_col_param is not greater than the number of subcolumns
        where(cv_col_param .gt. ncol) cv_col_param = ncol

      print *, 'cv_col_param(1:20) in gauss_fits=', cv_col_param(1:20)
      endif !if do_cv_param 
    end subroutine gauss_fits 
end module mod_gauss_fits

