import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def free_velocity_field(xy, t):
    x, w = xy
    vx = w
    vw = 0 
    return (vx, vw)


class multipart_simulation():
    def __init__(self,initial_condition,mass,t_span):
        
        self.fig,self.sub=plt.subplots(nrows=7,gridspec_kw={'height_ratios': [8, 1, 2,2,2,2,2]})
        self.mass=mass
        self.t_span=t_span
        self.dt=t_span[1]-t_span[0]
        #self.compute_simulation_free(np.array(initial_condition))
        self.compute_lyapunov_evol(np.array(initial_condition))
        self.stime = plt.Slider(self.sub[1], 'time', t_span[0], t_span[-1:], valinit=t_span[-1])
        self.stime.on_changed(self.updateplot)
        self.updateplot(t_span[-1])
        self.sub[2].plot(t_span,self.control[:,0])
        self.sub[3].plot(t_span,self.control[:,1])
        self.sub[4].plot(t_span,self.descent[:,0])
        self.sub[5].plot(t_span,self.descent[:,1])
        self.sub[6].plot(t_span,self.energy)
    def compute_simulation_free(self,initial_condition):
        self.sol=[]
        for xy0 in initial_condition:
            self.sol.append(odeint(free_velocity_field, xy0, self.t_span))
    def compute_lyapunov_evol(self,initial_condition):
        npart=len(initial_conditions)
        nt=len(self.t_span)
        self.sol=np.zeros((npart,nt,2))
        self.sol[:,0,:]=initial_conditions
        self.control=np.zeros((nt,2))
        self.descent=np.zeros((nt,2))
        self.energy=np.zeros(nt)
        for step in range(nt-1): 
            u_1=np.sum(-1*self.mass*self.sol[:,step,1]*np.cos(self.sol[:,step,0]))
            self.descent[step,0]=u_1
            u_1=np.sign(u_1)
            u_2=np.sum(-1*self.mass*self.sol[:,step,1]*np.sin(self.sol[:,step,0]))
            self.descent[step,1]=u_2
            u_2=np.sign(u_2)
            self.control[step,:]=u_1,u_2
            self.energy[step]=np.log(np.sum(0.5*self.mass*self.sol[:,step,1]*self.sol[:,step,1]))
            def F(xw,t):
                x,w=xw
                vx=w
                vw=u_1*np.cos(x)+u_2*np.sin(x)
                return (vx,vw)
            for k in range(npart):
                self.sol[k,step+1,:]=odeint(F,self.sol[k,step,:],[0,self.dt])[-1]



            
    def updateplot(self,t):
        step=int(t//(self.t_span[1]-self.t_span[0]))
        y=np.array([r[step,:] for r in self.sol])
        self.sub[0].clear()
        self.sub[0].scatter(np.mod(y[:, 0],2*np.pi), y[:, 1],s=1000*mass);
        self.sub[0].set_xlim(0, 2*np.pi)
        self.sub[0].set_ylim(np.min(self.sol[:,:,1]),np.max(self.sol[:,:,1]))
#initial_conditions =[(0,0.2),(0,0.3),(0,1)]
initial_conditions=[]
mass=[]
np.random.seed( 30 )
for i in range(25):
    for j in range(1):
        #initial_conditions.append((i*0.1,1+j*0.2))
        #initial_conditions.append((i*0.1,-1-j*0.2))
        #mass.append(1.1+np.sin(i))
        #mass.append(1.1+np.sin(i))
        initial_conditions.append(10*np.random.random(2)-5)
        mass.append(np.random.random())
mass=np.array(mass)
mass=mass/np.sum(mass)
print(mass)
t_span = np.arange(0, 400, 0.005)
sim=multipart_simulation(initial_conditions,mass,t_span)
plt.show()
#plt.axis('equal'); plt.xlabel('x'); plt.ylabel('y');
#plt.show()
