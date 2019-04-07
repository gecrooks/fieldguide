
set palette grey negative
#set terminal pdf enhanced monochrome solid size 3.5in, 2.1 font "/Library/Fonts/MyFonts/TrumpMediaevalLTStd-Roman.otf"
set terminal pdfcairo enhanced color transparent font "Trump Mediaeval LT Std,11" size 3.5in, 2.1in 
# background rgb 'gray'


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

BetaLogisticPDF(x,nu,lambda, A, G) = (1./Beta(A,G)) * abs(1./lambda)  * exp( -A * ((x-nu)/lambda)) * ( 1.+  exp( - ((x-nu)/lambda)) )**(-A-G)


set xtics 0,1
set ytics 0,0.5



PearsonVIIPDF(x, m) = ((1.+ (x**2) )**(-m)) / Beta(1./2, m-1./2)


#set xlabel 'ylabel' tc rgb 'white'
#set ylabel 'xlabel' tc rgb 'white'
#set border lc rgb 'white'
#set key tc rgb 'white'


StudentsTPDF(x, k) = ((1.+ (x**2)/k )**(-(k+1)/2.)) / (sqrt(k) * Beta(1./2, k/2.))

set style line 1 lt 1 lw 3 pt 3 lc rgb "gray"
unset label
set border 0
unset key; unset tics; unset border
set output "pdfCauchyNB.pdf"
set yrange[0:0.5]
plot[-4:4] \
	PearsonVIIPDF(x,1) lc rgb "white" lw 5

