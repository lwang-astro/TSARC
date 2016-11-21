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
For the propuse as we discussed above, we want to use a new variable s replacing the function of time \f$t\f$. 
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
Now we treat (t, Pt) as new coordinate and momemtum, H' as new Hamitonian, and use s as new integration variable.
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


Since \f$Pt = -H(t)\f$, \f$H'=H+Pt = T(\mathbf{P}) + U(\mathbf{Q}) = 0 \f$. Thus during integration, \f$T(\mathbf{P}) \approx -U(\mathbf{Q}) \f$. 
This requires \f$ f(T(\mathbf{P})) - f(-U(\mathbf{Q})) \approx 0 \f$. 
With Taylor expansion, we can obtain:

(13) \f$ f(T(\mathbf{P})) = f(-U(\mathbf{Q})) + \left[T(\mathbf{P}) + U(\mathbf{Q})\right] f'(-U(\mathbf{Q})) + O\left[T(\mathbf{P}) + U(\mathbf{Q})\right]^2 \f$

Thus \f$ g(\mathbf{Q},\mathbf{P}) \approx f'(-U(\mathbf{Q})) \f$


\subsubsection logH_sec Logarithmic Hamintonian method
Mikkola & Tanikawa (1999) suggests to use the function \f$ f(x) = \log{x} \f$ (Logarithmic Hamintonian method).

*/

