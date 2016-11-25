/*! \mainpage ARC Introduction
  
\section AR_sec Algorithmic Regularization (AR)

The algorithm used in this code is based on the literatures of <A HREF="http://adsabs.harvard.edu/abs/1999MNRAS.310..745M">Mikkola & Tanikawa (1999)</A> and <A HREF="http://adsabs.harvard.edu/abs/1999AJ....118.2532P">Preto & Tremaine (1999)</A>. The development of this code refers to the Chapter 2 (Mikkola) in book <A HREF="http://www.springer.com/us/book/9781402084300">The Cambridge N-body Lectures</A>. Here the basic idea of AR is described.

The numerical simulations of gravitational N-body systems dynamical evolutions are frequently used in astrophysics. 
However, due to the singularity of two-body gravitational potential when these two particles become infinite close lead to the difficulty in the highly accurately integration of two-body bounded system with very high eccentricity. 
To get high accuracy of integration when two particles are very close, the time step for integration should be reduced significantly.
This result in time consuming computation if number of particles is large.
On the other hand, for a long-term integration, the total energy of the systems may be systematiclly drifted due to the numerical accuracy of integrators. 
Thus the sympletic methods are suggested to be used since it can keep the energy conservation for long-term integration.

However, the sympletic methods are difficult to be applied for gravitational systems due to the required time step (integration step) shrinking when two particle get close.
Thus Mikkola & Tanikawa (1999) and Preto & Tremaine (1999) develop the special time transformation method based on the extended phase space Hamiltonian. 
The time \f$t\f$ become a general coordinate in Hamiltonian with corresponding general momentum \f$Pt\f$.
The integration of the equation of motion then depends on the new differential variable \f$ s\f$.
In this case, time and the motion of the system can be integrated with a fixed step size of s, which allow the usage of sympletic methods.

\subsection H_sec Hamiltonian in Extended Phase Space

Defining the general coordinates as \f$ \mathbf{q} = \{q_i\}, (i=1,n) \f$ with freedom of \f$n\f$ and corresponding general momentums ad \f$ \mathbf{p} \f$,  The Hamiltonian equations is:

(1) \f$ \frac{d \mathbf{q}}{d t} = \frac{\partial H}{\partial \mathbf{p}}\f$; \f$ \frac{d \mathbf{p}}{d t} = - \frac{\partial H}{\partial \mathbf{q}} \f$

Here the dt is used as a differetial varaible. 
For the propuse as we discussed above, we want to use a new variable \f$s\f$ replacing the function of time \f$t\f$. 
In this case, the time is treated as a new general coordinate. 
And the corresponding time momentum \f$Pt\f$ should be also added.

We extend the coordiantes to \f$ \mathbf{Q} = (t, \mathbf{q}) \f$ and the momentums to \f$ \mathbf{P} = (Pt, \mathbf{p})\f$  with total freedom of \f$2(n+1)\f$.

The new Hamiltonian \f$H'\f$ should also satisfy the Hamiltonian equations (1). 
Especially for \f$(t, Pt)\f$, we can get:

(2) \f$ \frac{d t}{d t} = \frac{\partial H'}{\partial Pt} = 1 \f$; \f$ \frac{d Pt}{d t} = - \frac{\partial H'}{\partial t} = - \frac{\partial H}{\partial t}\f$

From first equation of (2), we find the \f$H'\f$ linearly depend on \f$Pt\f$, thus \f$H'\f$ can be the form as \f$ H' = H + Pt \f$.
The second equation indicates that the time evolution of \f$Pt\f$ is equal to the negative energy change of the system.
Thus the value of \f$Pt\f$ at the time \f$t\f$ can be \f$-H(t)\f$.

We want to write Hamiltonian equations with new differetial variable \f$ ds\f$.
Defining \f$ g(\mathbf{Q},\mathbf{P}) = \frac{dt}{ds} \f$, we can rewrite (1) with \f$ds\f$ and extended coordinates \f$(\mathbf{Q}, \mathbf{P})\f$ as:

(3) \f$ \frac{d \mathbf{Q}}{d s} = g(\mathbf{Q},\mathbf{P}) \frac{\partial H'}{\partial \mathbf{P}} \f$;   \f$ \frac{d \mathbf{P}}{d s} = - g(\mathbf{Q},\mathbf{P}) \frac{\partial H'}{\partial \mathbf{Q}} \f$

However, we need to have the Hamiltonian equations the same form as original, thus we need to find another Hamiltonian \f$\Gamma(\mathbf{P},\mathbf{Q})\f$ that satisfy the Hamiltonian equations:

(4) \f$ \frac{d \mathbf{Q}}{d s} = \frac{\partial \Gamma}{\partial \mathbf{P}} \f$; \f$ \frac{d \mathbf{P}}{d s} = -\frac{\partial \Gamma}{\partial \mathbf{Q}} \f$

To find correct \f$\Gamma(\mathbf{P},\mathbf{Q})\f$, we go back to the Principle of least action which is used to derive the Lagrangian equations.
The relation between (standard) Hamiltonian \f$H(\mathbf{p},\mathbf{q},t)\f$ and Lagrangian \f$L(\mathbf{p},\mathbf{q},t)\f$ is 

(5) \f$ H(\mathbf{p},\mathbf{q},t) = \sum_{i=1}^n p_i \dot{q_i} - L(\mathbf{p},\mathbf{q},t) \f$

The Principle of least action require the action

(6) \f$ S = \int_{t_1}^{t_2} L(\mathbf{p},\mathbf{q},t) dt = \int_{t_1}^{t_2} \left[ \sum_{i=1}^n p_i \dot{q_i} - H(\mathbf{p},\mathbf{q},t) \right] dt \f$ 

should take the mimimum path, thus any function variation \f$ \delta S \f$ should makes \f$ S + \delta S\f$ increase.
Thus when \f$ \delta L(\mathbf{p},\mathbf{q},t) = 0 \f$, this condition is satisfied. This leads to the Lagrangian equations and also the Hamitonian equations.

Here the integration takes from \f$ t_1 \f$ to \f$ t_2 \f$ and the time is used as integration variable. 
Now we treat \f$(t, Pt)\f$ as new coordinate and momemtum, \f$H'\f$ as new Hamitonian, and use \f$s\f$ as new integration variable.
Then \f$S\f$ can be rewrited as:

(7) \f$ S = \int_{s_1}^{s_2} \left[ \sum_{i=1}^n p_i \frac{d q_i} {d s} + Pt \frac{d t}{d s} - (H(\mathbf{p},\mathbf{q},t) + Pt) \frac{d t}{d s} \right] ds = \int_{s_1}^{s_2} \left[ \sum_{i=1}^{n+1} P_i \frac{d Q_i}{d s} - H'(\mathbf{P},\mathbf{Q}) \frac{d t}{d s}\right] ds \f$

It is obvious that when

(8) \f$ \Gamma(\mathbf{P},\mathbf{Q}) = H'(\mathbf{P},\mathbf{Q}) \frac{d t}{d s} = g(\mathbf{Q},\mathbf{P}) (H(\mathbf{p},\mathbf{q},t) + Pt) \f$

The formula (7) become the same form as (6). Then with Principle of least action, the Hamiltonian equation (4) can be derived.
We call the \f$ \Gamma(\mathbf{P},\mathbf{Q}) \f$ is the Hamiltonian in the extended phase space \f$ (\mathbf{P},\mathbf{Q}) \f$

The Hamiltonian in extended phase space \f$\Gamma\f$ is also useful for analyzing the systems where Hamiltonian \f$H\f$ explicitly depends on time and is not conserved. 
Since time become a coordinate in \f$\Gamma\f$, \f$\frac{\partial \Gamma}{\partial s}\f$ is zero thus \f$ \Gamma\f$ become conserved quantity. 
The method dealing with closed system can be used with Hamiltonian in extended phase space.

\subsection T_sec Time transformation for Separable Hamiltonian

With the Hamiltonian in extended phase space, we can integrate the equation of motions with step \f$ ds \f$ by choosing a kind of \f$g(\mathbf{Q},\mathbf{P})\f$. If we choose a \f$g(\mathbf{Q},\mathbf{P})\f$ that makes the Hamiltonian \f$\Gamma(\mathbf{Q},\mathbf{P})\f$ separable for \f$P\f$ and \f$Q\f$:

(9) \f$ \Gamma(\mathbf{Q},\mathbf{P}) = a(\mathbf{P}) + b(\mathbf{Q}) \f$

Then explicit Leapfrog (SIA) integration method can be used. Preto & Tremaine (1999) suggests to use

(10) \f$ g(\mathbf{Q},\mathbf{P}) = \frac{f(T(\mathbf{P})) - f(-U(\mathbf{Q}))}{T(\mathbf{P}) + U(\mathbf{Q})} \f$

where \f$ T(\mathbf{P}) = T(\mathbf{p}) + Pt \f$ is the extended kinetic energy and \f$ U(\mathbf{Q}) = U(\mathbf{q},t) \f$ is the extended potential energy. 

The Hamiltonian becomes separable:

(11)  \f$ \Gamma = f(T(\mathbf{P})) - f(-U(\mathbf{Q})) \f$

Then the equation of motions are:

(12) \f$ \frac{d \mathbf{q} }{d s} = f'(T(\mathbf{p})+Pt) \frac{\partial T(\mathbf{p})}{\partial {\mathbf{p}}} \f$;
     \f$ \frac{d t }{d s} = f'(T(\mathbf{p})+Pt) \f$;
     \f$ \frac{d \mathbf{p} }{d s} = f'(-U(\mathbf{q},t)) \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{q}}} \f$;
     \f$ \frac{d Pt}{d s} = f'(-U(\mathbf{q},t)) \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{t}}} \f$;

where \f$ f'(x) = \frac{d f(x)}{d x} \f$.

Since \f$Pt = -H(t)\f$, \f$H'=H+Pt = T(\mathbf{P}) + U(\mathbf{Q}) = 0 \f$. Thus during integration, \f$T(\mathbf{P}) \approx -U(\mathbf{Q}) \f$. 
This requires \f$ f(T(\mathbf{P})) - f(-U(\mathbf{Q})) \approx 0 \f$. 
With Taylor expansion, we can obtain:

(13) \f$ f(T(\mathbf{P})) = f(-U(\mathbf{Q})) + \left[T(\mathbf{P}) + U(\mathbf{Q})\right] f'(-U(\mathbf{Q})) + O\left[T(\mathbf{P}) + U(\mathbf{Q})\right]^2 \f$

Thus 

(14) \f$ g(\mathbf{Q},\mathbf{P}) \approx f'(-U(\mathbf{Q})) \f$

\subsubsection logH_sec Logarithmic Hamintonian method (LogH)

Mikkola & Tanikawa (1999) suggests to use the function \f$ f(x) = \log{x} \f$ (Logarithmic Hamintonian method).
In this case, the time transformation based on (14) is:

(15) \f$ g(\mathbf{Q},\mathbf{P}) \approx \frac{1}{-U(\mathbf{Q})} \f$

Then the equation of motions can be written as:

(16) \f$ \frac{d \mathbf{q} }{d s} = \frac{1}{T(\mathbf{p})+Pt} \frac{\partial T(\mathbf{p})}{\partial {\mathbf{p}}} \f$;
     \f$ \frac{d t }{d s} = \frac{1}{T(\mathbf{p})+Pt} \f$;
     \f$ \frac{d \mathbf{p} }{d s} = \frac{1}{-U(\mathbf{q},t)} \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{q}}} \f$;
     \f$ \frac{d Pt}{d s} = \frac{1}{-U(\mathbf{q},t)} \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{t}}} \f$;

For the point mass systems with Newtonian gravity 

(17) \f$ T(\mathbf{p}) = \sum_{i=1}^{n} \frac{\mathbf{p_i}^2}{2m} \f$; \f$ U(\mathbf{q},t) = - \sum_{i<j,i=1,j=1}^{i\rightarrow n,j\rightarrow n} \frac{G m_i m_j}{|\mathbf{q_i}-\mathbf{q_j}|} \f$

where G is gravitational constant and \f$ m_i, m_j \f$ are masses of point-mass particles.

From (17) we see \f$ \frac{d Pt}{d s} = 0 \f$. 
This is only for the isolated system. If the system has external force from perturbers or external potential. The energy of system (\f$-Pt\f$) may not be conserved any more. Thus the energy change should be added into \f$Pt\f$ during the integration.

\subsubsection TTL_sec Time-Transformed Leapfrog (TTL)

The regularization methods where energy explicitly appear in the equation of motions cannot solve the few-body systems with large mass ratio (for example, planetary systems and super massive black hole with surrounding stars), because the energy is dominated by the massive bodies, and this introduce the systematic error during the integration. To solve this kind of issue, <A HREF="http://adsabs.harvard.edu/abs/2002CeMDA..84..343M">Mikkola & Aarseth (2002)</A> developed the so-called Time-Transformed Leapfrog (TTL) method. This method is also based on time transformation. The major difference compared with the LogH method is that the time transformation function also need to be integrated.

The time transformation (10) leads to the equations of motion (12) where time transformation \f$ f'(T(\mathbf{p})+Pt) \f$ and \f$ f'(-U(\mathbf{q},t))\f$ explicitly depend on kinetic energy, binding energy and potential. 
If we want to replace \f$ -U(\mathbf{q},t) \f$ to other quantity \f$ W(\mathbf{q})\f$, considering the requirement \f$ f(T(\mathbf{P})) - f(-U(\mathbf{Q})) \approx 0 \f$, we should also find another quantity \f$ w(\mathbf{p}) \f$ that allow \f$ f(w(\mathbf{p})) - f(W(\mathbf{q})) \approx 0 \f$. and 

(18) \f$ g(\mathbf{Q},\mathbf{P}) = \frac{f(w(\mathbf{p})) - f(W(\mathbf{q}))}{T(\mathbf{P}) + U(\mathbf{Q})} \approx f'(W(\mathbf{q})) \f$

Instead of finding the \f$ w(\mathbf{p}) \f$ for each kind of \f$ W(\mathbf{q})\f$, Mikkola & Aarseth (2002) suggest to use the differential equation 

(19) \f$  \frac{d W(\mathbf{q})}{d s} = \frac{\partial W(\mathbf{q})}{\partial \mathbf{q}} \cdot \frac{d \mathbf{q}} {d s} \f$

and integrate this equation to approximate \f$ w(\mathbf{p}) = \int \frac{d W(\mathbf{q})}{d s} d s\f$ simultaneously with integration of \f$ \frac{d \mathbf{p} }{d s} \f$.

However 

(20) \f$ \frac{d \mathbf{q}}{d s} = \frac{d \mathbf{q}}{d t} \frac{d t}{d s} = \frac {\mathbf{p}}{m} f'(W(\mathbf{q}))\f$

Thus \f$ \frac{d W(\mathbf{q})}{d s} \f$ explicitly depends on the momemtum. The integration in principle are not separatable. 
To solve this issue, Mikkole & Aarseth (2002) recommend to use averaged momemtums \f$ \langle \mathbf{p} \rangle \f$ (velocities) between previous and current step's during the Leapfrog integration, because the averaged values can represent the momemtums at the D (half) step when \f$\mathbf{q}\f$ is integrated.

Then if we take \f$ f(x) = \log{x}\f$ again, we have the equations of motion like:

(21) \f$ \frac{d \mathbf{q} }{d s} = \frac{1}{w} \frac{\partial T(\mathbf{p})}{\partial {\mathbf{p}}} \f$;
     \f$ \frac{d t }{d s} = \frac{1}{w} \f$;
     \f$ \frac{d \mathbf{p} }{d s} = \frac{1}{W(\mathbf{q})} \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{q}}} \f$;
     \f$ \frac{d w}{d s} = \frac{1}{W(\mathbf{q})} \frac{\partial W(\mathbf{q})}{\partial \mathbf{q}} \cdot \frac{\langle \mathbf{p} \rangle} {m} \f$;

This solution avoid use the energy (potential) as a time transformation dependence, thus with a suitable choice of \f$ W(\mathbf{q}) \f$, the high mass ratio systems can be integrated with high accuracy.

\section code_sec Implementation of ARC

We implememted AR method together with Chain (discussed below) for few-body systems by using C++ programming Language. 
The idea is make the integrator a C++ class thus can be easily used as a module for other codes. 
In this section we describe the details of the implementation.

\subsection chain_sec Particle Chain

If the bounded few-body systems are inside a big cluster enviroment, the average distance between these particles can be much smaller than the scale of cluster. 
Thus the round off error can be large if the positions of these particles are in the cluster center-of-mass frame.
To avoid this issue, <A HREF="http://adsabs.harvard.edu/abs/1993CeMDA..57..439M">Mikkola & Aarseth (1993)</A> suggested to use Chain method.

The idea is to connect all particles in one chain and using relative position and velocity for integration.
Firstly, one particle is selected as a starting point of the chain, then the nearest particle is selected as the next chain member, the relative position \f$ X \f$ and velocity \f$ V \f$ between these neighbors are calculated and stored.
After that, we found the nearest particle to this second member from the remaining particles and calculate relative positions and velocites and do this iterately until all particles are connected in this chain.
The relative positions and velocites can be described by absolute positions and velocities in a ordered chain as:

(22) \f$ \mathbf{X}_i = \mathbf{q}_{i+1} - \mathbf{q}_i \f$; \f$ \mathbf{V}_i = \mathbf{v}_{i+1} - \mathbf{v}_i \f$

The integration is done with these relative quantities to reduce round off error. The equations of motion can be written as 

(23) \f$ \frac{d \mathbf{X}_i}{d t} = \mathbf{V}_i \f$; \f$ \frac{d \mathbf{V}_i}{d t} = \mathbf{A}_{i+1} - \mathbf{A}_i \f$

where \f$ \mathbf{A}_i \f$ is the acceleration of particle \f$ i\f$.

When the particles are moved, the nearest neighbor of each particle may become different, thus the update of chain order should be performed with a suitable time interval. 

\subsection leap_sec Leapfrog Integrator

By combining the AR algorithm and Chain scheme, we can construct a Leapfrog integrator of equations of motion for $N$-body systems like:
- D mode:

(24)      \f$ \Delta t = \Delta s / (\alpha (T(\mathbf{p}) + Pt) + \beta w + \gamma) \f$;
          \f$ t += \Delta t \f$;
          \f$ \mathbf{X}_i += \Delta t \mathbf{V}_i \f$ 

- K mode:

(25)      \f$ \delta t = \Delta s / (\alpha U(\mathbf{q},t) + \beta W(\mathbf{q}) + \gamma) \f$;
          \f$ \mathbf{V}_i += \delta t (\mathbf{A}_{i+1} - \mathbf{A}_{i}) \f$;
          \f$ Pt += \delta t \sum_i (-m_i \langle \mathbf{v}_i \rangle \cdot f_{ext,i}) \f$;
          \f$ w += \delta t \sum_i \frac{\partial W}{\partial \mathbf{q}_i} \cdot \langle \mathbf{v}_i \rangle \f$ 

where \f$ f_{ext,i} \f$ is the external force from outside the system (e.g., perturber force or tidal force) of each particle \f$ i\f$, and \langle \mathbf{v}_i \rangle is obtained by averaging the velocities of the initial and the final \f$ \mathbf{v}_i \f$ of this K mode step. \f$ \alpha, \beta, \gamma \f$ are the coefficients representing the weights of the LogH, TTL and non-time-transformation modes separately. For example, if \f$ \alpha=0\f$, then no LogH will be performed, and if \f$ \alpha =1, \beta=0, \gamma=0 \f$ it is LogH ARC.

The initial value of \f$ Pt \f$ should be the initial binding energy of the system \f$ U(\mathbf{q},t) - T(\mathbf{p}) \f$. 
If the system is isolated, \f$ Pt \f$ is constant.
The initial value of \f$ w\f$ is set to initial \f$ W(\mathbf{q}) \f$.

The Leapfrog step start with half-step D and then loop full-step K-D-K and stop with half-step D:

(26) \f$ D(s/2)K(s)D(s)....K(s)D(s/2) \f$

This provide a second order integrator of ARC. Trying this integrator for a two-body bounded system can result in an energy and eccentricity conserved kepler orbit. Only the time phase can have cumulative error after long-term integration.

\subsection extrapolation_sec Extrapolation Integrator

The Leapfrog integrator only has second order accuracy, which is not enough for many applications. One can reduce the step size of integration to obtain higher accuracy. 
However, as energy is always conserved for two-body motions, we don't have good checker to indicate whether the integration is accurate enough. 
A better and more efficient way is to extrapolate the integration results to obtain high order accuracy.
The idea of extrapolation integration is well summarized in <A HREF="http://link.springer.com/book/10.1007%2F978-0-387-21738-3"> Stoer & Bulirsch</A>. 
Here the basic algorithm is shown.

First, if we integrate the equations of motion with Leapfrog integrator by step \f$ s\f$. 
we get the first result with a certain accuracy. Now we keep the total step constant but divide the integration into several sub-steps with equal sizes by \f$ n \f$, we can obtain higher accuracy of the integration. When we use a sequence of dividers \f$ (n_1, n_2, n_3 ...)\f$ (\f$ n_{i+1}>n_i\f$) and do the integration with each of them, we can obtain a series of results with increasing accuracy. Then we can extrapolate these results to get an even higher accurate one.

There are two major methods of extrapolation: polynomial and rational.
Both methods can be described as recursive functions:

- Polynomial: 

(27) \f$ T_{i,k} = T_{i,k-1} + \frac{T_{i,k-1} - T_{i-1,k-1}}{( h_{i-k} / {h_i} )^2 -1} \f$, \f$ 1 \le k \le i \le m \f$ 

- Rational:

(28) \f$ T_{i,k} = T_{i,k-1} + \frac{T_{i,k-1} - T_{i-1,k-1}}{\left[ \frac{h_{i-k}}{h_i} \right]^2 \left[ 1 - \frac{T_{i,k-1} - T_{i-1,k-1}}{T_{i,k-1}- T_{i-1,k-2}} \right]-1} \f$, \f$ 1 \le k \le i \le m \f$ 

Here \f$ i\f$ indicate the integration with sub-step size \f$ h_i = s/n_i\f$, and \f$ k \f$ indicate the extrapolation order.
The \f$ T_{i,0} \f$ are results of Leapfrog integrations, and for each order \f$ i\f$, the \f$ T_{i,i} \f$ is final extrapolation result we want. 
The \f$ T_{i,i} \f$ can be obtained by calculating \f$ T_{i,k} \f$ from \f$ k=1 \f$ to \f$ k=i \f$ using the recursive functions.

One benefit of these recursive functions is that a higher order extrapolation \f$ T_{i+1,i+1} \f$ can be established based on current existing \f$ T_{i,k}, k=0\sim i \f$  with a new higher order integration result \f$ T_{i+1,0} \f$. 
Then it is easy to estimate the error by comparing \f$ T_{i+1,i+1} \f$ and \f$ T_{i,i} \f$ to determine whether another higher order result is necessary. 
For example, in ARC integration, we can check the time or position phase error and energy error to determine how many orders we need to integrate and extrapolate due to the accuracy requirment.

The sequences of dividers \f$ n_i \f$ have several choices for different applications:
- <A HREF="https://en.wikipedia.org/wiki/Romberg's_method">Romberg</A>: (1, 2, 4, 8 ...)
- <A HREF="http://link.springer.com/article/10.1007%2FBF02165234">Bulirsch & Stoer</A> (BS): (1, 2, 3, 4, 6, 8 ...)
- <A HREF="http://link.springer.com/article/10.1007%2FBF01385634">Hairer</A> (4k): (2, 6, 10, 14 ...)
- Harmonic:  (1, 2, 3, 4 ...)

Different seuqnces and recursive functions can be combined together for extrapolation integration. 
We implement all sequences shown above. Later we discuss the special application of some sequences.

\subsection dense_sec Dense Output for Time Synchronization


\subsection step_sec Integration Step Control

\subsection perf_sec Performance Analysis
          
*/
/*
Consider the systems with positions \f$\mathbf{r}\f$, velocity \f$\mathbf{v}\f$ and acceleration \f$\mathbf{F}(\mathbf{r})\f$ and the time transformation function \f$ \frac{d t}{d s} = \frac{1}{\Omega(\mathbf(r))} \f$. The equation of motions can be described as:

(18) \f$ \frac{d \mathbf{r}} {d s} = \frac{\mathbf{v}}{\Omega(\mathbf(r))} \f$; \f$ \frac{d t}{d s} = \frac{1}{\Omega(\mathbf(r))}\f$ ; \f$ \frac{d \mathbf{v}}{d s} = \frac{\mathbf{F}(\mathbf{r})}{\Omega(\mathbf(r))} \f$ 

As discussed in 
*/
