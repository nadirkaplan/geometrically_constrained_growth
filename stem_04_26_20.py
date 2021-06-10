#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 12:10:40 2020

@author: nadir
"""


#from fenics import *
from precipitation_general import *
#import matplotlib.pyplot as plt
#plt.ion()
#import sys, math
#import numpy as np
#import time
from datetime import date

set_log_level(30)

kk=1
delta=0.1
mm=sqrt(0.1)     
sigmamin=Domain(kk).min
sigmamax=Domain(kk).max
chi=Bend(0.02).chi

        
#time parameters 
T =int(round(45)) #int(round(25))            # final time
dt=Constant(1e-7) #1e-6, 1e-7 are stable.
dtarray=np.zeros([1]) #this will initiate a numpy array with the first element being the initial value of dt.
dtarray[0]=float(dt)   #the index nn in the main loop will also correspond to the indices of dtarray.
dtmax=1e-3 #maximum step size

class Rate: #parameters of the growth rate
    def __init__(self):
    
        #coral
        self.alpha1=Constant(1)
        self.alpha2=Constant(1)
        self.alpha3=Constant(0.5)
        self.eta1=Constant(1.6)
        self.eta21=Constant(-0.3)
        self.eta22=Constant(0.5)
        self.eta3=Constant(0.)
        self.eta41=Constant(0.)
        self.eta42=Constant(0.)
        self.eta5=Constant(0.)
        self.eta6=Constant(0.)
        self.lambda1=Constant(1.)
        
        '''
        self.alpha1=Constant(1)
        self.alpha2=Constant(1)
        self.alpha3=Constant(0.5)
        self.eta1=Constant(1.)
        self.eta21=Constant(0.2)
        self.eta22=Constant(-0.1)
        self.eta3=Constant(1.)
        self.eta41=Constant(1.)
        self.eta42=Constant(1.)
        self.eta5=Constant(0.)
        self.eta6=Constant(0.)
        self.lambda1=Constant(0.4)

        '''


# tfi=0.01  #below tfi, time step is dt1
# tff=T    #between tfi and tff, time step is dt2
# dt1= Constant(1e-5)
# dt2= Constant(1e-3)   #Constant(1e-3); 5e-3 timestep is unstable.
# dt = dt1 # initial time step size; this will be updated in the time loop

# ntime1=int(round(tfi/dt1))
# ntime2=int(round((tff-tfi)/dt2))
# ntotal=ntime1+ntime2
# num_steps=ntotal
# dtarray=np.zeros([ntotal])
# dtarray[0:ntime1]= dt1
# dtarray[ntime1:ntotal]= dt2

nframe=100;
#deltat_frame=T/nframe;
#tframe=float(dt1);
toutput=np.zeros([nframe])
toutput[:]=[T/nframe*kk for kk in range(nframe)]
toutput[0]=float(dtarray[0])
toutput[-1]=float(T)

#create mesh and define function space
mesh=IntervalMesh(5000, float(sigmamin), float(sigmamax))
#mesh.smooth(10)
xmesh=mesh.coordinates()


#xarraytemp=[1e-8+(0.25-1e-8)/6000*jj for jj in range (0, 6001)] record for post processing

tau=1.0 #1.5
xmesh[:]=geomseq(xmesh, float(sigmamin), float(sigmamax), tau)


P1 = FiniteElement('CG', interval, 2)
P2 = FiniteElement('DG', interval, 2)

Neqs=int(6)
NeqsAux=int(9)
element1 = MixedElement([P1, P1, P2, P2, P2, P2]) #MixedElement([P1]*Neqs)
element2 = MixedElement([P2, P2, P2, P2, P2, P2, P2, P2, P2])

V = FunctionSpace(mesh, element1)
Q = FunctionSpace(mesh, element2)

# Define test functions

v1, v2, v3, w1, w2, w3 = TestFunctions(V)
y1, y2, y3, p1, p2, p3, h1, h2, h3 = TestFunctions(Q)

u = Function(V)
un = Function(V)
un_1 = Function(V)
un_2 = Function(V)
un_3 = Function(V)
un_4 = Function(V)
du = TrialFunction(V)

q = Function(Q)
qn = Function(Q)
dq = TrialFunction(Q)

#initial value expressions in String format

gg0='( 2*pi+2*delta*mm*mm*pi*cos(2*kk*pi*x[0])+pow(delta*kk*mm*mm, 2.)*pi*pow(sin(2*kk*pi*x[0]), 2.) )'


kn10temp=('( 1+mm*mm*(-0.5+delta*(-1 +kk*kk)*cos(2*kk*pi*x[0]))+(pow(mm, 4.)*(3+delta*delta*(4-6*kk*kk)'
          +'-4*delta*(-1+kk*kk)*cos(2*kk*pi*x[0])+2*delta*delta*(2-5*kk*kk)*cos(4*kk*pi*x[0])))/8. )')

kn10='('+kn10temp+')'

kg0='( -mm+pow(mm, 3.)*(0.5 - delta*(-1 + kk*kk)*cos(2*kk*pi*x[0])) )'

taug0='0.'

kn20='0.'

H0='('+kn10+'+'+kn20+')/2.'


YY10='( cos(2*pi*x[0])*(1+delta*mm*mm*cos(2*kk*pi*x[0])) )'
YY20='(-1-delta*mm*mm*cos(2*kk*pi*x[0]))*sin(2*pi*x[0])'
YY30='0.'

tau10temp=('( -sin(2*pi*x[0])-delta*kk*mm*mm*cos(2*pi*x[0])*sin(2*kk*pi*x[0])'
           +'+(delta*delta*kk*pow(mm, 4.)*sin(2*kk*pi*x[0])*(2*cos(2*pi*x[0])*cos(2*kk*pi*x[0])'
           +'+kk*sin(2*pi*x[0])*sin(2*kk*pi*x[0])))/2. )')
tau20temp=('( -cos(2*pi*x[0])+delta*kk*mm*mm*sin(2*pi*x[0])*sin(2.*kk*pi*x[0])'
           +'+(delta*delta*kk*pow(mm, 4.)*sin(2*kk*pi*x[0])*(-2*cos(2*kk*pi*x[0])*sin(2*pi*x[0])'
           +'+kk*cos(2*pi*x[0])*sin(2*kk*pi*x[0])))/2. )')
tau30temp='0.'

tau10=tau10temp
tau20=tau20temp
tau30=tau30temp

NN10temp=('(-cos(2*pi*x[0])+mm*mm*(cos(2*pi*x[0])/2+delta*kk*sin(2*pi*x[0])*sin(2*kk*pi*x[0]))'+
      '+(pow(mm, 4.)*(cos(2*pi*x[0])*(-3+2*pow(delta*kk, 2.)-2*pow(delta*kk, 2.)*cos(4*kk*pi*x[0]))'+
      '-4*delta*kk*(1+2*delta*cos(2*kk*pi*x[0]))*sin(2*pi*x[0])*sin(2*kk*pi*x[0])))/8.)')
NN20temp=('(sin(2*pi*x[0])+mm*mm*(-sin(2*pi*x[0])/2+delta*kk*cos(2*pi*x[0])*sin(2*kk*pi*x[0]))'+
      '+(pow(mm, 4.)*((3 - 2*pow(delta*kk, 2.)+2*pow(delta*kk, 2.)*cos(4*kk*pi*x[0]))*sin(2*pi*x[0])'+
      '-4*delta*kk*cos(2*pi*x[0])*(1+2*delta*cos(2*kk*pi*x[0]))*sin(2*kk*pi*x[0])))/8.)')
NN30temp=('(mm-pow(mm, 3.)/2.)')

NN10=NN10temp
NN20=NN20temp
NN30=NN30temp


#Do not change these!
kG0=kn10+'*'+kn20+'-'+taug0+'*'+taug0 

U01=('-alpha1*'+kg0+'+alpha2*'+kg0+'*'+kg0+'+alpha3*'+kg0+'*'+kg0+'*'+kg0)
U02=('+eta1*'+H0+'*'+H0+'+eta21*'+kn10+'*'+kn20
     +'+eta22*'+taug0+'*'+taug0+'+eta3*'+kg0+'*'+H0+'*'+H0+'-eta41*'+kg0+'*'+kn10+'*'+kn20
     +'+eta42*'+kg0+'*'+taug0+'*'+taug0+'+eta5*'+taug0+'*'+H0+'+eta6*'+kg0+'*'+taug0+'*'+H0)
U0=U01+U02


#Define initial values
#(kg0, U0, kn20, gg0, kn10, taug0)
u0=Expression((kg0, U0, kn20, gg0, kn10, taug0), 
              degree=1, delta=delta, kk=kk, mm=mm, alpha1=Rate().alpha1, alpha2=Rate().alpha2, 
              alpha3=Rate().alpha3, eta1=Rate().eta1, eta21=Rate().eta21, eta22=Rate().eta22,
              eta3=Rate().eta3, eta41=Rate().eta41, eta42=Rate().eta42, eta5=Rate().eta5, eta6=Rate().eta6, pi=pi)

q0=Expression((YY10, YY20, YY30, tau10, tau20, tau30, NN10, NN20, NN30), 
              degree=1, delta=delta, kk=kk, pi=pi, mm=mm)

un=interpolate(u0, V)
qn=interpolate(q0, Q)

kg, U, kn2, gg, kn1, taug = split(u);
kgn, Un, kn2n, ggn, kn1n, taugn = split(un);
kgn_1, Un_1, kn2n_1, ggn_1, kn1n_1, taugn_1 = split(un_1);
kgn_2, Un_2, kn2n_2, ggn_2, kn1n_2, taugn_2 = split(un_2);
kgn_3, Un_3, kn2n_3, ggn_3, kn1n_3, taugn_3 = split(un_3);
kgn_4, Un_4, kn2n_4, ggn_4, kn1n_4, taugn_4 = split(un_4);

#plot(Un)


YY1, YY2, YY3, tau1, tau2, tau3, NN1, NN2, NN3 = split(q);
YY1n, YY2n, YY3n, tau1n, tau2n, tau3n, NN1n, NN2n, NN3n = split(qn);

#Variables for flux boundary conditions

g10=Constant(0.)   #dkg/dsigma at sigma=0
g20=Constant(0.)   #dU/dsigma  at sigma=0
g30=Constant(0.)   #dkn2/dsigma at sigma=0

g11=Constant(0.)   #dkg/dsigma at sigma=1/2
g31=Constant(0.)   #dkn2/dsigma at sigma=1/2

boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)

bx0 = BoundaryX0(kk)
bx0.mark(boundary_markers, 0)

bx1 = BoundaryX1(kk)
bx1.mark(boundary_markers, 1)
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

param=Rate()
lambdau=param.lambda1
gammab=Bend().gamma
zetab=Bend().zeta


#BDF coefficients:

ccn=Constant(1.)
ccn_1=Constant(1.)
ccn_2=Constant(1.)
ccn_3=Constant(1.) #These values are only for initiation, they'll be updated accordingly in the main cell.
ccn_4=Constant(1.)
ccr= Constant(1.) #weight of the entire RHS.



#new weak formulation for the main computation at O(delta t^3) accuracy
#we will need to worry about nonzero flux boundary conditions later


F1=inner(lambdau*grad(kg)/gg, grad(v1))*dx-gg*U*v1*dx+gg*v1*Utemp(kg, kn1, kn2, taug, param, H)*dx

F2=((gg*kg+ccn*ggn*kgn+ccn_1*ggn_1*kgn_1+ccn_2*ggn_2*kgn_2+ccn_3*ggn_3*kgn_3+ccn_4*ggn_4*kgn_4)/dt*v2*dx
    +inner(ccr*Q2(U, gg), grad(v2))*dx-ccr*f2(U, gg, kg, kn1, kn2, taug, kappaG)*v2*dx)

#F3=((gg*kn2+ccn*ggn*kn2n+ccn_1*ggn_1*kn2n_1+ccn_2*ggn_2*kn2n_2+ccn_3*ggn_3*kn2n_3+ccn_4*ggn_4*kn2n_4)/dt*v3*dx
#    +inner(ccr*Q3(gg, kn2), grad(v3))*dx-ccr*f3(U, gg, kg, kn1, kn2, qbend)*v3*dx)

F3=((H(kn1, kn2)+ccn*H(kn1n, kn2n)+ccn_1*H(kn1n_1, kn2n_1)+ccn_2*H(kn1n_2, kn2n_2)+ccn_3*H(kn1n_3, kn2n_3)+ccn_4*H(kn1n_4, kn2n_4))/dt*v3*dx
    -ccr*f3(chi, U, kn1, kn2, H)*v3*dx)

F4=((gg+ccn*ggn+ccn_1*ggn_1+ccn_2*ggn_2+ccn_3*ggn_3+ccn_4*ggn_4)/dt*w1*dx-ccr*f4(U, gg, kg)*w1*dx)

F5=((gg*kn1+ccn*ggn*kn1n+ccn_1*ggn_1*kn1n_1+ccn_2*ggn_2*kn1n_2+ccn_3*ggn_3*kn1n_3+ccn_4*ggn_4*kn1n_4)/dt*w2*dx
    -ccr*f5(U, gg, kg, kn1, kn2, taug)*w2*dx)

F6=((gg*taug+ccn*ggn*taugn+ccn_1*ggn_1*taugn_1+ccn_2*ggn_2*taugn_2+ccn_3*ggn_3*taugn_3+ccn_4*ggn_4*taugn_4)/dt*w3*dx
    -ccr*f6(U, gg, kg, kn1, kn2, taug)*w3*dx)


#Auxiliary equations:

nn1=NN2*tau3-NN3*tau2
nn2=NN3*tau1-NN1*tau3
nn3=NN1*tau2-NN2*tau1

F7=( (YY1-YY1n)/dt-nn1*U )*y1*dx
F8=( (YY2-YY2n)/dt-nn2*U )*y2*dx
F9=( (YY3-YY3n)/dt-nn3*U )*y3*dx

F10=( gg*(tau1-tau1n)/dt-gg*U*taug*NN1-U.dx(0)*nn1 )*p1*dx
F11=( gg*(tau2-tau2n)/dt-gg*U*taug*NN2-U.dx(0)*nn2 )*p2*dx
F12=( gg*(tau3-tau3n)/dt-gg*U*taug*NN3-U.dx(0)*nn3 )*p3*dx

F13=( (NN1-NN1n)/dt+U*kn2*nn1+U*taug*tau1 )*h1*dx
F14=( (NN2-NN2n)/dt+U*kn2*nn2+U*taug*tau2 )*h2*dx
F15=( (NN3-NN3n)/dt+U*kn2*nn3+U*taug*tau3 )*h3*dx


F=F1+F2+F3+F4+F5+F6
Faux=F7+F8+F9+F10+F11+F12+F13+F14+F15

#main

nmesh=len(mesh.coordinates())

#MODIFY ALL OF THESE!
#ntime=round(num_steps/100)
uarray=np.zeros([nframe, Neqs, nmesh])
qarray=np.zeros([nframe, NeqsAux, nmesh])
timearray=np.zeros([nframe])

shape='lantern'


today = date.today()
d3 = today.strftime("%y_%m_%d")
filenamemesh=shape+'_mesh_kk_'+str(kk)+'_'+d3+'.txt'
filenameX=shape+'_Xdata_kk_'+str(kk)+'_'+d3+'.txt'
filenameY=shape+'_Ydata_kk_'+str(kk)+'_'+d3+'.txt'
filenameZ=shape+'_Zdata_kk_'+str(kk)+'_'+d3+'.txt'
filenametime=shape+'_timearray_kk_'+str(kk)+'_'+d3+'.txt'

filenamekg=shape+'_kappag_kk_'+str(kk)+'_'+d3+'.txt'
filenameU=shape+'_U_kk_'+str(kk)+'_'+d3+'.txt'
filenamekN2=shape+'_kappaN2_kk_'+str(kk)+'_'+d3+'.txt'
filenamegg=shape+'_gg_kk_'+str(kk)+'_'+d3+'.txt'
filenamekN1=shape+'_kappaN1_kk_'+str(kk)+'_'+d3+'.txt'
filenametaug=shape+'_taug_kk_'+str(kk)+'_'+d3+'.txt'

tic = time.perf_counter()

dv = TrialFunction(V)
   
#initial guess computation

L=F
a=-derivative(L, u, du)
du=Function(V)
dkg, dU, dkn2, dgg, dkn1, dtaug = split(du);
u.assign(un)

BDFparameters(0, ccn, ccn_1, ccn_2, ccn_3, ccn_4, ccr)

tol = 1.0E-7
iter = 0
maxiter = 25
eps = 1.0
# u_k must have right boundary conditions here

while eps > tol and iter < maxiter:
    iter += 1
    print(iter, 'iteration')
    
    #bc_dU0 = DirichletBC(V.sub(1), -u(float(sigmamin))[1], bx0)
    #bc_dU = DirichletBC(V.sub(1), -u(float(sigmamax))[1], bx1)
    #bc_dtaug= DirichletBC(V.sub(5), -u(float(sigmamax))[5], bx1)
    A, b = assemble_system(a, L) #, [bc_dU0, bc_dU]) #, bc_dtaug)
    solve(A, du.vector(), b)    
    eps = np.linalg.norm(du.vector()[:], ord=np.Inf)
    print('Norm:', eps)
    u.assign(u+du)


solve(Faux == 0, q) 


#end of initial guess computation



#bc_kn2 = DirichletBC(V.sub(2), Constant(0.), bx0) # boundaryleft
#bc_taug= DirichletBC(V.sub(5), Constant(0.), bx1) #boundaryright
#bcs = [bc_taug] #trivially satisfied


# Time-stepping
t = dtarray[0]
kk= 0
nn=0
Nit=2; Nitmax=4; Nitmin=1;

#for nn=0:
    
#assign(u.sub(1), un.sub(1))

qn.assign(q)    
BDFparameters(nn, ccn, ccn_1, ccn_2, ccn_3, ccn_4, ccr)
BDFvariables(nn, u, un, un_1, un_2, un_3, un_4)
J = derivative(F, u, dv)  
#Jq= derivative(Faux, q, dq)



while  t<T: # and nn<500:
    
    
    #print(nn)
    #print('dtarray:', dtarray)

    # Solve variational problem for time step
    #solve(F == 0, u, bcs) #, solver_parameters={"newton_solver":{"relative_tolerance":1e-6}}
    
    problem = NonlinearVariationalProblem(F, u, J=J) # bcs, J)
    solver  = NonlinearVariationalSolver(problem)
    
    
    
    prm = solver.parameters
    prm['nonlinear_solver']='newton'
    prm["newton_solver"]["absolute_tolerance"]= 1e-7   #1e-9  5921 secs
    prm["newton_solver"]["relative_tolerance"] = 1e-6 #1e-7
    prm["newton_solver"]["maximum_iterations"] = 100
    prm["newton_solver"]["error_on_nonconvergence"]=False
    
    #prm['newton_solver']['linear_solver'] = 'mumps'
    #info(solver.parameters, True)
    
    
    (Nit, conv)=solver.solve()
    if conv==False:
        print(conv)
        print('time:', t,'; index:', nn)
        break
    
    
    #auxiliary equations:
    solve(Faux == 0, q) #, J=Jq)
    
    #solve(F == 0, u, [bc_u1, DirichletBC(V.sub(1), u3, bx1), bc_u3], solver_parameters={"newton_solver":{"relative_tolerance":1e-6}})

    # Save solution to file (VTK)
    #_u2, _u3 = u.split()
    #vtkfile_u2 << (_u2, t)
    #vtkfile_u3 << (_u3, t)

    # Update previous solution
    #un.assign(u)
    
       # Update previous solution
       
    #print('time:', float(t), 'time step:', float(dt), 'iterations:', Nit)
    
    if Nit<Nitmax:
        
               
        if t<=toutput[kk]<(t+float(dt)):
        
            print('index:', nn, '; iterations:', Nit, 'time:', float(t), 'time step:', float(dt))
        
            uarray[kk, :, :]= [ [u(ii)[jj] for ii in mesh.coordinates()] for jj in range(Neqs) ]
            qarray[kk, :, :]= [ [q(ii)[jj] for ii in mesh.coordinates()] for jj in range(NeqsAux) ]
                 
            timearray[kk]=t
            kk+=1
        
        qn.assign(q)
        
        if dtarray[nn]<(dtmax/2.):
            dt.assign(Constant(2.*dtarray[nn]))
        elif dtmax/2.<dtarray[nn]<dtmax:
            dt.assign(Constant(dtmax))
                   
        nn+=1
        BDFparameters(nn, ccn, ccn_1, ccn_2, ccn_3, ccn_4, ccr)
        BDFvariables(nn, u, un, un_1, un_2, un_3, un_4)
        
            
        J = derivative(F, u, dv)     
        #Jq= derivative(Faux, q, dq)
        
        # Update current time
        t += float(dt)
        dtarray=np.append(dtarray, float(dt))
        
        #print('dtarray:', dtarray)
        
        
    else:
        dt.assign(Constant(dtarray[nn]/2.))
        dtarray[nn]=dtarray[nn]/2.
        print('time:', float(t), 'time step:', float(dt))
        #dtarray[-1]=float(dt)
    
    #if n%100==0: 
   
    



        
toc = time.perf_counter()    

print(f"Total runtime is {toc - tic:0.4f} seconds.")



X1array=qarray[0:kk, 0, :];  X2array=qarray[0:kk, 1, :]; X3array=qarray[0:kk, 2, :];
kappag_array=uarray[0:kk, 0, :]; U_array=uarray[0:kk, 1, :]; kN2_array=uarray[0:kk, 2, :];
gg_array=uarray[0:kk, 3, :]; kN1_array=uarray[0:kk, 4, :]; taug_array=uarray[0:kk, 5, :];



np.savetxt(filenamekg, kappag_array.transpose(),  fmt='%.6e', delimiter='\t')
np.savetxt(filenameU, U_array.transpose(),  fmt='%.6e', delimiter='\t')
np.savetxt(filenamekN2, kN2_array.transpose(),  fmt='%.6e', delimiter='\t')
np.savetxt(filenamegg, gg_array.transpose(),  fmt='%.6e', delimiter='\t')
np.savetxt(filenamekN1, kN1_array.transpose(),  fmt='%.6e', delimiter='\t')
np.savetxt(filenametaug, taug_array.transpose(),  fmt='%.6e', delimiter='\t')

np.savetxt(filenamemesh, mesh.coordinates(),  fmt='%.6e')
np.savetxt(filenameX, X1array.transpose(),  fmt='%.6e', delimiter='\t')
np.savetxt(filenameY, X2array.transpose(),  fmt='%.6e', delimiter='\t')
np.savetxt(filenameZ, X3array.transpose(),  fmt='%.6e', delimiter='\t')
np.savetxt(filenametime, timearray,  fmt='%.6e', delimiter='\t')


