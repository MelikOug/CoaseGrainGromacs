using BlackBoxOptim
using Optim
using NumericalIntegration
using Plots
using CSV	
using DataFrames
using Roots
using NLsolve
import LinearAlgebra: norm
#All units have been put in units of k_BT at 300K; at some point if temperature studies important should do this explicitly

#Read in the parameters from args

#The order shall be  beta(in kBT) EpsA(unitless) EpsB(unitless) rsolventA(nm) rsolventB(nm) muA(D) muB(D) alphaA(D^2/eV) alphaB(D^2/eV) IPA(eV) IPB(eV) dipolepolymer(D) polarisabilitypolymer(D^2/eV) IPpolymer (eV) rPolymer(nm) R(nm) PhiAve(unitless) Name(string)

#Solvent interaction terms ⍺ the ratio of volumes between solvent A and B and the β the interaction
β=parse(Float64, ARGS[1])

#Eps of solvents used in bruggeman
εA=parse(Float64, ARGS[2])
εB=parse(Float64, ARGS[3])

#Average radius of solvent molecules. (nm)
rA=parse(Float64, ARGS[4])
rB=parse(Float64, ARGS[5])

#Dipole Moments of Solvents (D)
μA=parse(Float64, ARGS[6])
μB=parse(Float64, ARGS[7])

#Polarisability of Solvents (D^2/eV)

⍺A=parse(Float64, ARGS[8])
⍺B=parse(Float64, ARGS[9])

#IPs for the solvents (eV)

IPA=parse(Float64, ARGS[10])
IPB=parse(Float64, ARGS[11])

#Number density of solvent molecules at 300K
ρA=3/(4*pi*rA^3)
ρB=3/(4*pi*rB^3)
println("RhoA is ",ρA," whilst RhoB is ", ρB)

⍺v=ρA/ρB

#Polymer Attributes
μP=parse(Float64, ARGS[12])
⍺P=parse(Float64, ARGS[13])
IPP=parse(Float64, ARGS[14])
rP=parse(Float64, ARGS[15])

R=parse(Float64,ARGS[16])

PhiAve=parse(Float64, ARGS[17])

minrS=min(rA,rB)
resstep=0.01 #in nm

NaT=ρA * PhiAve * 4*pi*((R+resstep)^3-(rP+minrS)^3)/3
NbT=ρB * (1-PhiAve) * 4*pi*((resstep+R)^3-(rP+minrS)^3)/3

Name=ARGS[18]
println("Sucessfully read in input:")
println("Volumetric ratio of solvents is: ", ⍺v)
println("Cost of mixing of solvents parameter is: ", β)
println("Properties of the polymer: dipole moment=",μP, "D, polarisability=", ⍺P, " IP=",IPP," and rp=",rP)
println("Properties of the solvent A: dipole moment=",μA, "D, polarisability=", ⍺A, " IP=",IPA," and rp=",rA)
println("Properties of the solvent B: dipole moment=",μB, "D, polarisability=", ⍺B, " IP=",IPB," and rp=",rB)
println("The system as specified as a volume of ",  4*pi*(R^3-rP^3)/3," nm^3 with ", NaT," solvent A molecules and ",NbT," solvent B molecules in it.")
# Set the resolution of all the calculations
#Keesom Constant (1/(24*pi^2*eps0^2)) in units of eV^2 nm^6/D^4
KeesomConst=0.000000259709403
KT=0.0257

#Calculate ⍺ Ua-Ub
UA=0.0 .- (36 .* KeesomConst .* ⍺P .* ⍺A  .* IPP .* IPA ./ ( IPP .+ IPA ) ) .- (KeesomConst .* μP .^ 2 .* μA .^ 2 ./ KT ) .- ( 1.5 .* KeesomConst .* μP .^2 .* ⍺A )
UB=0.0 .- (36 .* KeesomConst .* IPP .* IPB ./ ( IPP .+ IPB ) .* ⍺P .* ⍺B ) .- (KeesomConst .* μP .^ 2 .* μB .^ 2 ./ KT ) .- ( 1.5 .* KeesomConst .* μP .^2 .* ⍺B )

#Now we can calculate the constant in alpha*Ua-Ub propto PotentialConst/r^6
PotentialConst=((⍺v .* UA) .- UB) ./ KT

println("Calculated the potential term components to be UA=",UA," UB=",UB," thus constant to be: ", PotentialConst)
#Quantities used in rcutoff energy induced calulatione
k0= 18.09512802   #1/eps0 in \frac{eV nm}{e^2}
Delta= (ρA .* μA .- ρB .* μB) .* 0.02081943
B=ρB .*  μB .* 0.02081943 

function SecondDerivCheck(ɸ,alpha,beta)
	return  alpha* (- alpha^3 * (1 + 2* beta) * (ɸ - 1)^2 + ɸ * (2 * beta * ( ɸ -1 ) + ɸ ) - alpha^2* (ɸ -1) * (1 + (-3 + 4*beta)*ɸ) + alpha * (4* beta* (ɸ - 1)^2 + (2 - 3* ɸ)* ɸ) + 4*(alpha -1)* alpha * (alpha* (ɸ-1) - ɸ)*(ɸ-1) * atanh((alpha - (1 + alpha) * ɸ)/(alpha + ɸ - alpha * ɸ))) / ((ɸ - 1) * (alpha + ɸ * (1 - alpha ))^4)
end

function EquationOfMotionForVolumeFraction(ɸ,r,Const,Power,SolventRelativeDensity,SolventInteractionEnergy,Lambda, Na,Nb)
	return abs(Const/r^Power + SolventRelativeDensity*log(SolventRelativeDensity*ɸ[1]/((SolventRelativeDensity-1)*ɸ[1]+1))-log((1-ɸ[1])/((SolventRelativeDensity-1)*ɸ[1]+1)) + SolventRelativeDensity * SolventInteractionEnergy * (1-2*ɸ[1] -(SolventRelativeDensity-1)*ɸ[1]^2) / ((SolventRelativeDensity-1)*ɸ[1] + 1)^2 + Lambda )
end



function EffectiveEpislonEquation(EpsEff, PolymerRadius, BoxRadius,EnergytoOptimiseTo)
    return abs.(EnergytoOptimiseTo .- (PolymerRadius .^ -2 .- BoxRadius .^ -2) .* (EpsEff[1]-1).* (2 .* EpsEff[1]+1) ./ EpsEff[1] .^ 3 ./ 2 )
end
 
function LambdaConstraint(Radius,VolFrac,PolymerRadius,RadiusBox,BoxTotalVolFrac)
    return  integrate(Radius,r.*r.*VolFrac)[]-BoxTotalVolFrac*(RadiusBox^3-PolymerRadius^3)/3
end


function LambdaFunctionToOpt(Lambda,Const,Power,SolventRelativeDensity,SolventInteractionEnergy,PolymerRadius,RadiusBox,BoxTotalVolFrac,NaT,NbT,RhoA,RhoB)
    r=range(PolymerRadius,RadiusBox,step=resstep)
    Φ=zeros(length(r),1)
    OutputNa=zeros(length(r),1)
    OutputNb=zeros(length(r),1)
    for i in eachindex(r)
	    
	if (i >2)
		Na=OutputNa[i-1] - RhoA* 4*pi*Φ[i-1]*(r[i-1].^3-r[i-2].^3)/3#RhoA*integrate(r[1:(i-1)],Φ[1:(i-1)] .* 4 .* pi .* r[1:(i-1)] .^ 2)[]
		Nb=OutputNb[i-1] - RhoB* 4*pi*(1-Φ[i-1])*(r[i-1].^3-r[i-2].^3)/3#RhoB*integrate(r[1:(i-1)],(1 .- Φ[1:(i-1)]) .* 4 .* pi .* r[1:(i-1)] .^ 2)[]
	#	println("r",r[i-1]," Na is ", Na," NB is ",Nb)
		OutputNa[i]=Na
		OutputNb[i]=Nb
		if(Na<0)
			Na=0.0
                end
                if(Nb<0)
                        Nb=0.0
                end 
	else
		Na=NaT
		Nb=NbT
		OutputNa[i]=Na
		OutputNb[i]=Nb
	end
	res=bboptimize(p->EquationOfMotionForVolumeFraction(p,r[i],Const,Power,SolventRelativeDensity,SolventInteractionEnergy,Lambda[],Na,Nb),[0.5];NumDimensions=1,SearchRange=(0.0,1.0),TraceMode=:silent,Method=:de_rand_1_bin_radiuslimited) #:adaptive_de_rand_1_bin_radiuslimited)
	Φ[i]=best_candidate(res)[]

    end
    err=(integrate(r,r.*r.*Φ)[]-BoxTotalVolFrac*(RadiusBox^3-PolymerRadius^3)/3)
    outputdata=[r Φ]
    outputfile1=string("lambdamodel",Name, ".dat")
    CSV.write(outputfile1,DataFrame(outputdata,:auto),delim=' ')
    return err
end

		


#lresult0=optimize(λ->GuessLambdaFunctionToOpt(λ,PotentialConst,6,⍺v,β,rP,R,PhiAve), [0.0], BFGS())
#λ=bisection(λ->LambdaFunctionToOpt(λ,PotentialConst,6,⍺v,β,rP,R,PhiAve),-5.0,19);
#lresult=optimize(λ->LambdaFunctionToOpt(λ,PotentialConst,6,⍺v,β,rP,R,PhiAve), [Optim.minimum(lresult0)], BFGS())
minrS=min(rA,rB)
λ=find_zero(λ->LambdaFunctionToOpt(λ,PotentialConst,6,⍺v,β,rP+minrS,R,PhiAve,NaT,NbT,ρA,ρB), -1.5,Order0(),  verbose=true)
println("find_zero algorithm found λ to be ",λ)

#println("Final result for lambda", Optim.minimum(lresult))
#lresult=bboptimize(λ->LambdaFunctionToOpt(λ,PotentialConst,6,⍺v,β,rP,R,PhiAve);NumDimensions=1,Method=:dxnes, MaxTime = 100)
#println("Final result for lambda using bboptimize", best_candidate(lresult)[])
#λ=best_candidate(lresult)[]
#λ=find_zero(λ->LambdaFunctionToOpt(λ,PotentialConst,6,⍺v,β,rP,R,PhiAve), 0.1, Order0(),verbose=true)
roundedstartr=round((rP+minrS) / resstep) * resstep
r=range(roundedstartr,R,step=resstep)
Φ=zeros(length(r),1)
x=zeros(length(r),1)
ε=zeros(length(r),1)
ERR=zeros(length(r),1)
MinCheck=zeros(length(r),1)
energy=zeros(length(r),1)
outputNA=zeros(length(r),1)
outputNB=zeros(length(r),1)
for i in eachindex(r)
	if (i >2)
	     Na=outputNA[i-2] - ρA* 4*pi*Φ[i-1]* (r[i-1].^3-r[i-2].^3)/3   #integrate(r[1:(i-1)],Φ[1:(i-1)] .* 4 .* pi .* r[1:(i-1)] .^ 2)[]
	     Nb=outputNB[i-2] - ρB* 4*pi*(1-Φ[i-1])* (r[i].^3-r[i-2].^3)/3 #ρB*integrate(r[1:(i-1)],(1 .- Φ[1:(i-1)]) .* 4 .* pi .* r[1:(i-1)] .^ 2)[]
             #if(Na<0)
	     #	     Na=0.0
	     #end
	     #if(Nb<0)
	     #	     Nb=0.0
	     #end
        else
	     Na=NaT
             Nb=NbT
	     outputNA[i]=Na
	     outputNB[i]=Nb
        end
    outputNA[i]=Na
    outputNB[i]=Nb
    res=bboptimize(p->EquationOfMotionForVolumeFraction(p,r[i],PotentialConst,6,⍺v,β,λ,Na,Nb),[0.5];NumDimensions=1,SearchRange=(0.0,1.0),TraceMode=:silent,Method=:de_rand_1_bin_radiuslimited) #:adaptive_de_rand_1_bin_radiuslimited)

	#	    res=bboptimize(p->EquationOfMotionForVolumeFraction(p,r[i],PotentialConst,6,⍺v,β,λ),[Φ[i-1]];NumDimensions=1,SearchRange=(0.0,1.0),TraceMode=:silent,Method=:xnes)
    Φ[i]=best_candidate(res)[]
    ERR[i]=best_fitness(res)[]
    MinCheck[i]=SecondDerivCheck(Φ[i],⍺v,β)
    x[i]=ρA*Φ[i]/(Φ[i]*(ρA-ρB)+ρB)
end

#Φculm= integrate(r, Φ .* 4 .* pi .* r .* r) ./(r.^3 .- (rP+minrS-0.00001) .^     3)

#Calculate the volume fraction of the sphere at r from it at shell at r.
Φculm= cumul_integrate(r, Φ .* 4 .* pi .* r .* r) ./cumul_integrate(r,4 .* r .* r.* pi)
for i in eachindex(r)
    #Bruggeman Solution using volume fraction to calculate permittivity at each r
    ε[i]=((3*Φculm[i]-1)*εA+(2-3*Φculm[i])*εB)/4+sqrt(((3*Φculm[i]-1)*εA+(2-3*Φculm[i])*εB)^2+8*εA*εB)/4
end


#We now have the eps the model predicts but we to calculate an effective epsilon

#We start by calculating the total induced energy in the system of the solvent on the polymer

OnPolymerUA=0.0 .- (36 .* KeesomConst .* ⍺P .* ⍺A  .* IPP .* IPA ./ ( IPP .+ IPA ) ) .- (KeesomConst .* μP .^ 2 .* μA .^ 2 ./ KT ) .- ( 1.5 .* KeesomConst .* μA .^2 .* ⍺P )
OnPolymerUB=0.0 .- (36 .* KeesomConst .* IPP .* IPB ./ ( IPP .+ IPB ) .* ⍺P .* ⍺B ) .- (KeesomConst .* μP .^ 2 .* μB .^ 2 ./ KT ) .- ( 1.5 .* KeesomConst .* μB .^2 .* ⍺P )


Delta=(ρA .* μA .^2 .- ρB .* μB .^2 ) ./ 3.14159 .^ 2 .* 0.02081943 .^ 2 ./ KT .^2 ./ 3
B=ρB .* μB .^ 2 ./ 3.14159 .^ 2 * 0.02081943 .^2 ./ KT .^ 2 ./ 3


Eindscreenedtotal=integrate(r, ((ε .- 1) .* (2 .* ε .+ 1 )) ./ (r .^ 3 .* ε .^ 3)   )

IndPercentagescreenedall=  0.0 .- ⍺P .* 0.02081943 .^ 2 .* k0 .^2 ./ 18 .* (cumul_integrate(r, ((ε.- 1) .* (2 .* ε .+ 1 ) ) ./ (r .^ 3 .* ε .^ 2)  )) .^ 2 ./ 18 ./ Eindscreenedtotal .* 100




#Now we need to solve this equal to the energy induced by a single effective medium of epseff


εeffscreened = optimize(p -> EffectiveEpislonEquation(p, rP+minrS, R,Eindscreenedtotal[1]),2.0,80.0,[2.1], Fminbox(LBFGS()),Optim.Options(show_trace=false, allow_f_increases=true, iterations=2000, f_tol=1e-15, x_tol=1e-15) ; autodiff = :forward).minimizer[]

println("Epsilon from full model equating with screening= ",εeffscreened)

println("The error in the governing equatons for eps with screening is: ",EffectiveEpislonEquation(εeffscreened, rP, R,Eindscreenedtotal[1]))
#This finds the Rcutoff initially then the model predicted Epsilon

indexscreened=findmin(abs.(ε .- εeffscreened))[2]

rcutoffscreened=r[indexscreened]
volfracoutscreened=Φculm[indexscreened]
xfracoutscreened=volfracoutscreened*ρA/(volfracoutscreened*(ρA-ρB)+ρB)
εoutputscreened=ε[indexscreened]

println("rcutoff with screening ", rcutoffscreened)
println("molefrac with screening ", xfracoutscreened)
println("molefrac whole system ", PhiAve*ρA/(PhiAve*(ρA-ρB)+ρB))
println("Preferential Solvation ",xfracoutscreened-PhiAve*ρA/(PhiAve*(ρA-ρB)+ρB))
println("volfrac with screening ", volfracoutscreened)
println("epsilon at rcutoff with screening ", εoutputscreened)

energysystemtot=integrate(r,4 .* pi .* r.^2 .* energy)
println("Energy of Total System",energysystemtot)
#When the difference drops below thermal noise (1KT) then new solvent no longer has an impact

outputdata=[r x Φ Φculm ε outputNA outputNB ERR MinCheck energy]
outputparams=[PotentialConst λ rcutoffscreened εoutputscreened energysystemtot xfracoutscreened-PhiAve*ρA/(PhiAve*(ρA-ρB)+ρB)] #These are in units of KT, KT, nm and unitless

outputfile1=string("model",Name, ".dat")
outputfile2=string("model",Name, "_computedparams.dat")

CSV.write(outputfile1,DataFrame(outputdata,:auto),delim=' ')
CSV.write(outputfile2,DataFrame(outputparams,:auto))

finalerr=(integrate(r,r.*r.*Φ)[]-PhiAve*(last(r)^3-first(r)^3)/3)
println("Final Error in Constraint is ", finalerr)

 
