

#set terminal pdf enhanced monochrome solid size 3.5in, 2.1 font "/Library/Fonts/MyFonts/TrumpMediaevalLTStd-Roman.otf"
set terminal pdfcairo enhanced monochrome solid font "Trump Mediaeval LT Std,11" size 3.5in, 2.1in 


set samples 1001

# theta 
# alpha A
# beta B
# gamma G

Beta(A,G) = ((gamma(A)*gamma(G))/ gamma(A+G) )

AmorosoPDF(x, a, theta, alpha, beta) =(sgn(theta) * x <a) ? 0 : (1./gamma(alpha) ) * abs(beta/theta) * ((x-a)/theta)**(alpha*beta-1) * exp(-((x-a)/theta)**(beta) )

GammaPDF(x,theta,alpha) = AmorosoPDF(x,0.,theta,alpha,1.)

GammaExpPDF(x,nu,lambda,alpha) = (1./gamma(alpha)) * abs(1./lambda) * exp(-alpha * ((x-nu)/lambda) - exp(-(x-nu)/lambda) )



UnitGammaPDF(x, a, s, A, B) = (abs(B)/( gamma(A)*s )) * ((x-a)/s)**(B-1.) * ( - B * log ( (x-a)/s) )**(A-1.)


GenBetaPDF(x, a, s, A, G,B) = (gamma(A+G)/ (gamma(A)*gamma(G)) ) * abs(B/s) * ((x-a)/s)**(A*B-1.) * (1.-((x-a)/s)**B)**(G-1.)

BetaPrimePDF(x, a, s, A, G,B) = (gamma(A+G)/ (gamma(A)*gamma(G)) ) * abs(B/s) * ((x-a)/s)**(A*B-1.) * (1.+((x-a)/s)**B)**(-A-G)

BetaPDF(x,a,s,A,G) = GenBetaPDF(x,a,s,A,G,1.)

NormalPDF(x,a,s) = exp( - ( (x-a)**2 / (2.*s**2) )  ) / (sqrt (2.*pi* s**2) )

LogNormalPDF(x,a,V,B) = (V/(x-a)) * exp(- (B * log((x-a)/V))**2 /2) * abs(B) / (sqrt(2 * pi * V**2)) 



BetaExpPDF(x,nu,lambda, A, G) = (1./Beta(A,G)) * abs(1./lambda)  * exp( -A * ((x-nu)/lambda)) * ( 1.-  exp( - ((x-nu)/lambda)) )**(G-1.)



set xtics 0,1
set ytics 0,0.5



set output "pdfGammaPDF.pdf"
set border 3
set xtics nomirror
set ytics nomirror
# set title "Gamma(1/{/Symbol a}, {/Symbol a})  (unit variance)"
set label "{/Symbol a}=1" at first 0.1, first 0.95
set label "{/Symbol a}=2" at first 0.3, first 0.8
#set label "{/Symbol a}=3" at first 0.35, first 0.85
set label "{/Symbol a}=4" at first 0.4, first 0.9
#set label "{/Symbol a}=5" at first 0.45, first 0.95
set label "{/Symbol a}=6" at first 0.5, first 1.0
set label "{/Symbol a}=8" at first 0.55, first 1.1
set nokey
plot [0:3] \
	GammaPDF(x,1.,1.) lt 1, \
	GammaPDF(x,1./2,2.) lt 1, \
	GammaPDF(x,1./4,4.) lt 1, \
	GammaPDF(x,1./6,6.) lt 1, \
	GammaPDF(x,1./8,8.) lt 1

unset label
set border 3
set output "pdfAmorosoBetaPDF.pdf"
set nokey
#set title "Amoroso(x,0,1,2,{/Symbol b})"
set label "{/Symbol b}=4" at first 1.4, first 1.4
set label "{/Symbol b}=3, Wilson-Hilferty" at first 1.5, first 1.1
set label "{/Symbol b}=2, scaled chi" at first 1.8, first 0.6
set label "{/Symbol b}=1, gamma" at first 2.3, first 0.4
plot [0:3] \
	AmorosoPDF(x,0,1.,2.,4.) lt 1, \
	AmorosoPDF(x,0,1.,2.,3.) lt 1, \
	AmorosoPDF(x,0,1.,2.,2.) lt 1, \
	AmorosoPDF(x,0,1.,2.,1.) lt 1

unset label
set border 3
set output "pdfAmorosoBetaNegPDF.pdf"
set nokey
#set title "Amoroso(x,0,1,2,{/Symbol b}), negative {/Symbol b}"
set label "{/Symbol b}=-1\n inverse\n gamma " at first 0.06, first 1.6
set label "{/Symbol b}=-2\n scaled\n inverse-chi" at first 0.25, first 1.90
set label "{/Symbol b}=-3" at first 0.55, first 2.
plot [0:2] \
	AmorosoPDF(x,0,1.,2.,-1.) lt 1, \
	AmorosoPDF(x,0,1.,2.,-2.) lt 1, \
	AmorosoPDF(x,0,1.,2.,-3.) lt 1


unset label
set border 3
#set logscale y
set output "pdfStretchedExpPDF.pdf"
set nokey
#set title "Amoroso(x|0,1,1,{/Symbol b}), stretched exponential"
set label "{/Symbol b}=1" at first 0.1, first 1.0
set label "{/Symbol b}=2" at first 0.35, first 0.85
set label "{/Symbol b}=3" at first 0.5, first 1.0
set label "{/Symbol b}=4" at first 0.55, first 1.2
set label "{/Symbol b}=5" at first 0.6, first 1.4
plot [0:3] \
	AmorosoPDF(x,0,1.,1,1.) lt 1, \
	AmorosoPDF(x,0,1.,1,2.) lt 1, \
	AmorosoPDF(x,0,1.,1,3.) lt 1, \
	AmorosoPDF(x,0,1.,1,4.) lt 1, \
	AmorosoPDF(x,0,1.,1,5.) lt 1

unset label
set border 3
set output "pdfAmorosoBeta2PDF.pdf"
set nokey
#set title "Amoroso(x|0,1,{/Symbol a}, 2)"
set label "{/Symbol a}=1/2, half-normal" at first 0.3, first 1.1
set label "{/Symbol a}=1, Rayleigh" at first 0.7, first 0.9
set label "{/Symbol a}=3/2, Maxwell" at first 1.3, first 0.8
plot [0:3] \
	AmorosoPDF(x, 0, 1., 0.5, 2.) lt 1, \
	AmorosoPDF(x, 0, 1., 1.0, 2.) lt 1, \
	AmorosoPDF(x, 0, 1., 1.5, 2.) lt 1

unset label
set border 3
set output "pdfChiSqr.pdf"
set nokey
#set title "ChiSqr(x|k)"
set yrange [0:0.5]
set label "k=1" at first 0.6, first 0.45
set label "k=2" at first 0.25, first 0.30
set label "k=3" at first 0.3, first 0.25
set label "k=4" at first 0.5, first 0.16
set label "k=5" at first 0.8, first 0.05
plot [0:8] \
	AmorosoPDF(x, 0, 2, 0.5, 1.) lt 1, \
	AmorosoPDF(x, 0, 2, 1.0, 1.) lt 1, \
	AmorosoPDF(x, 0, 2, 1.5, 1.) lt 1, \
	AmorosoPDF(x, 0, 2, 2.0, 1.) lt 1, \
	AmorosoPDF(x, 0, 2, 2.5, 1.) lt 1

unset label
set border 3
set output "pdfChiP.pdf"
set nokey
set label "inverse chi, k=1" at first 0.6, first 1.4
set label "inverse chi square, k=1" at first 0.6, first 1.4
set label "chi, k=1" at first 0.6, first 1.4
set label "chi square, k=1" at first 0.6, first 1.4
set yrange [0:1.5]
plot [0:3] \
	AmorosoPDF(x, 0, 1./2, 1., -1.) lt 1, \
	AmorosoPDF(x, 0, 1./sqrt(2.), 1., -2.) lt 1, \
	AmorosoPDF(x, 0, 2., 1., 1.) lt 1,  \
	AmorosoPDF(x, 0, sqrt(2.), 1., 2.)  lt 1

unset label
set border 3
set output "pdfEVD.pdf"
set yrange[0:1]
set xtics -3,1
#set title "Extreme value distributions"
set label "standard Gumbel" at first -0.6, first 0.90
set label "reversed Weibull, {/Symbol b}=2" at first -2.9, first 0.8
set label "Frechet, {/Symbol b}=-2" at first 1.2, first 0.8
plot[-3:3] \
	GammaExpPDF(x,0,1,1) lt 1, \
	AmorosoPDF(x,0,1,1,-2) lt 1, \
	AmorosoPDF(x,0,-1,1,2) lt 1


unset label
set border 3
set output "pdfGammaExp.pdf"
set yrange[0:1]
set label "{/Symbol a}=1" at first 0.5, first 0.4
set label "{/Symbol a}=2" at first -0.1, first 0.5
set label "{/Symbol a}=3" at first -0.5, first 0.6
set label "{/Symbol a}=4" at first -0.9, first 0.75
set label "{/Symbol a}=5" at first -1.2, first 0.87
plot[-3:3] \
	GammaExpPDF(x,0,1,1) lt 1, \
	GammaExpPDF(x,0,1,2) lt 1, \
	GammaExpPDF(x,0,1,3) lt 1, \
	GammaExpPDF(x,0,1,4) lt 1, \
	GammaExpPDF(x,0,1,5) lt 1


unset label
set border 3
set output "pdfGumbel.pdf"
set yrange[0:0.5]
plot[-4:8] \
	GammaExpPDF(x,0,1,1) lt 1


PearsonVIIPDF(x, m) = ((1.+ (x**2) )**(-m)) / Beta(1./2, m-1./2)

unset label
set border 3
set output "pdfCauchy.pdf"
set yrange[0:0.5]
plot[-4:4] \
	PearsonVIIPDF(x,1) lt 1

StudentsTPDF(x, k) = ((1.+ (x**2)/k )**(-(k+1)/2.)) / (sqrt(k) * Beta(1./2, k/2.))

unset label
set border 3
set output "pdfStudentsT.pdf"
set yrange[0:0.5]
plot[-4:4] \
	StudentsTPDF(x,1.) lt 1, \
	StudentsTPDF(x,2.) lt 1, \
	StudentsTPDF(x,3.) lt 1, \
	NormalPDF(x,0,1) lt 1
	

unset label
set border 11
set output "pdfBeta.pdf"
set yrange[0:3]
plot[0:1] \
	BetaPDF(x,0.,1.,1.,1.) lt 1, \
	BetaPDF(x,0.,1.,2.,2.) lt 1, \
	BetaPDF(x,0.,1.,8.,8.) lt 1


unset label
set border 3
set output "pdfStdExp.pdf"
set yrange[0:1.1]
plot[0:4] exp(-x)


unset label
set border 3
set output "pdfNormal.pdf"
set yrange[0:1.0]
set xtics -4,2
set label "{/Symbol s}=2" at first 1.0, first 0.6
set label "{/Symbol s}=1" at first 1.5, first 0.3
set label "{/Symbol s}=0.5" at first 3.0, first 0.15
plot[-6:6] NormalPDF(x,0,1) lt 1, NormalPDF(x,0,0.5) lt 1, NormalPDF(x,0,2.) lt 1


unset label
set border 3
set output "pdfLogNormal.pdf"
set yrange[0:1.8]
set xtics 0,0.5
set label "{/Symbol b}=1" at first 0.08, first 0.75
set label "{/Symbol b}=2" at first 0.35, first 0.9
set label "{/Symbol b}=4" at first 0.50, first 1.2
plot[0.01:3]  \
	LogNormalPDF(x,0,1,1) lt 1 , \
	LogNormalPDF(x,0,1,2) lt 1, \
	LogNormalPDF(x,0,1,4) lt 1



unset label
set border 11
set output "pdfUnitGamma.pdf"
set yrange[0:3.0]
set xtics 0,0.5
set label "{/Symbol a}=1.5, {/Symbol b}=1" at first 0.08, first 2
set label "{/Symbol a}=2, {/Symbol b}=2" at first 0.25, first 1.6
set label "{/Symbol a}=5, {/Symbol b}=8" at first 0.33, first 2.1
plot[0.01:1]  UnitGammaPDF(x,0,1,1.5,1) lt 1, UnitGammaPDF(x,0,1,2,2) lt 1, UnitGammaPDF(x,0,1,5,8) lt 1


unset label
set border 3
set output "pdfUnitGammaNegBeta.pdf"
set yrange[0:1.0]
set xtics 0,0.5
set label "{/Symbol a}=1.5, {/Symbol b}=-1" at first 1.30, first 0.225
set label "{/Symbol a}=2, {/Symbol b}=-1" at first 1.40, first 0.36
set label "{/Symbol a}=5, {/Symbol b}=-8" at first 2.0, first 0.7
plot[1:4]  UnitGammaPDF(x,0,1,1.5,-1) lt 1, UnitGammaPDF(x,0,1,2,-1) lt 1, UnitGammaPDF(x,0,1,5,-8) lt 1





unset label
set border 11
set output "pdfPowerFn.pdf"
PowerFnPDF(x,b) = abs(b) * x ** (b-1.)
set yrange[0:4.0]
set xtics 0,0.2
set label "{/Symbol b}" at first 1.02, first 1.25
set label "1.5" at first 1.0075, first 1.5
set label "2.0" at first 1.0075, first 2.0
set label "2.5" at first 1.0075, first 2.5
set label "3.0" at first 1.0075, first 3.0
set label "3.5" at first 1.0075, first 3.5
set label "4.0" at first 1.0075, first 4.0
plot[0:1]   PowerFnPDF(x,1.5) lt 1,  PowerFnPDF(x,2) lt 1,  PowerFnPDF(x,2.5) lt 1, PowerFnPDF(x,3) lt 1,  PowerFnPDF(x,4) lt 1, PowerFnPDF(x,3.5) lt 1

unset label
set border 11
set output "pdfPowerFn2.pdf"
set yrange[0:4.0]
set xtics 0,0.2
set label "{/Symbol b}" at first 1.02, first 1.25
set label "0.2" at first 1.0075, first 0.2
set label "0.4" at first 1.0075, first 0.4
set label "0.6" at first 1.0075, first 0.6
set label "0.8" at first 1.0075, first 0.8
plot[0:1]    PowerFnPDF(x,0.2) lt 1, PowerFnPDF(x,0.4) lt 1, PowerFnPDF(x,0.6) lt 1,  PowerFnPDF(x,0.8) lt 1

unset label
set border 3
set output "pdfPareto.pdf"
set yrange[0:4.0]
set xtics 0,0.5
plot[1:3]   PowerFnPDF(x,-1.) lt 1,  PowerFnPDF(x,-2.) lt 1,  PowerFnPDF(x,-4.) lt 1



unset label
set border 3
set output "pdfLaplace.pdf"
set yrange[0:0.6]
set xtics -3.,1.
plot[-3:3]  exp(- abs(x ) )/2.





unset label
set border 3
set output "pdfBetaExp.pdf"
set yrange[0:1.]
set xtics -3.,1.
plot[0:4]   BetaExpPDF(x,0., 1., 2., 2.0) lt 1, BetaExpPDF(x,0., 1., 2., 4.0) lt 1, BetaExpPDF(x,0., 1., 2., 8.0) lt 1



unset label
set border 3
set output "pdfExpExp.pdf"
set yrange[0:1]
set xtics -3.,1.
plot[0:4]  BetaExpPDF(x,0., 1., 1., 2.)



unset label
set border 3
set output "pdfSinhNK.pdf"
set yrange[0:1]
set xtics -3.,1.
plot[0:4]  BetaExpPDF(x,0., 0.5, 0.25, 0.5) lt 1, BetaExpPDF(x,0., 1., 0.5, 0.5) lt 1


LogLogisticPDF(x,a,s,B) = BetaPrimePDF(x, a, s, 1., 1.,B)
unset label
set border 3
set output "pdfLogLogistic.pdf"
set yrange[0:2.1]
plot[0:2]  LogLogisticPDF(x,0,1.,1) lt 1,LogLogisticPDF(x,0,1.,2) lt 1,LogLogisticPDF(x,0,1.,4) lt 1, LogLogisticPDF(x,0,1.,8) lt 1, LogLogisticPDF(x,0,1.,0.5) lt 1



