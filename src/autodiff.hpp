#ifndef AUTODIFF_HPP
#define AUTODIFF_HPP

#include <array>  
#include <cmath> // for sin and cos

namespace ASC_ode
{

  template <size_t N, typename T = double>
  class Variable 
  {
    private:
      T m_val;
    public:
      Variable (T v) : m_val(v) {}
      T value() const { return m_val; }
  };

  template <typename T = double>
  auto derivative (T v, size_t /*index*/) { return T(0); } 


  template <size_t N, typename T = double>
  class AutoDiff
  {
  private:
    T m_val;
    std::array<T, N> m_deriv;
  public: 
    AutoDiff () : m_val(0), m_deriv{} {}
    AutoDiff (T v) : m_val(v), m_deriv{} 
    {
      for (size_t i = 0; i < N; i++)
        m_deriv[i] = derivative(v, i);
    }
    
    template <size_t I>
    AutoDiff (Variable<I, T> var) : m_val(var.value()), m_deriv{} 
    {
      m_deriv[I] = 1.0;
    }

    T value() const { return m_val; }
    std::array<T, N>& deriv() { return m_deriv; }
    const std::array<T, N>& deriv() const { return m_deriv; }
  };


  template <size_t N, typename T = double>
  auto derivative (AutoDiff<N, T> v, size_t index) 
  {
    return v.deriv()[index];
  }



  template <size_t N, typename T>
  std::ostream & operator<< (std::ostream& os, const AutoDiff<N, T>& ad)
  {
    os << "Value: " << ad.value() << ", Deriv: [";
    for (size_t i = 0; i < N; i++)
    {
      os << ad.deriv()[i];
      if (i < N - 1) os << ", ";
    }
    os << "]";
    return os;
  }

  template <size_t N, typename T = double>
  AutoDiff<N, T> operator+ (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
  {
     AutoDiff<N, T> result(a.value() + b.value());
     for (size_t i = 0; i < N; i++)
        result.deriv()[i] = a.deriv()[i] + b.deriv()[i];
       return result;
   }

   template <size_t N, typename T = double>
   auto operator+ (T a, const AutoDiff<N, T>& b) { return AutoDiff<N, T>(a) + b; }


   template <size_t N, typename T = double>
   AutoDiff<N, T> operator* (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
   {
       AutoDiff<N, T> result(a.value() * b.value());
       for (size_t i = 0; i < N; i++)
          result.deriv()[i] = a.deriv()[i] * b.value() + a.value() * b.deriv()[i];
       return result;
   }

   using std::sin;
   using std::cos;

   template <size_t N, typename T = double>
   AutoDiff<N, T> sin(const AutoDiff<N, T> &a)
   {
       AutoDiff<N, T> result(sin(a.value()));
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = cos(a.value()) * a.deriv()[i];
       return result;
   }

   // Exercise 18.4.
   // 1. Add additional useful operators for the AutoDiff class

   // New. - operator (a-b)
   template <size_t N, typename T=double>
   AutoDiff<N, T> operator-(const AutoDiff<N, T>&a, const AutoDiff<N, T>&b)
   {
       AutoDiff<N, T> result(a.value() - b.value()); // save the original values in result
       // now we use the math derivative expression (a'-b') and we save it in result
       for (size_t i = 0; i < N; i++) 
           result.deriv()[i] = a.deriv()[i] - b.deriv()[i];
       return result;
   }

   // New. negative value operator (-b). It is not included in the above because 0-x are 2
   // variables but when you write -a it is only readed as 1 var, not like 0-a, so we implement
   // its derivative here in a separate function with only 1 input var
   template <size_t N, typename T=double>
   AutoDiff<N, T> operator-(const AutoDiff<N, T> &a)
   {
       AutoDiff<N, T> result(-a.value());
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = -a.deriv()[i];
       return result;
   }

   // New. Division Operator. (a/b)'=((a'*b-a*b')/b^2)
   template <size_t N, typename T=double>
   AutoDiff<N, T> operator/(const AutoDiff<N, T> &a,  const AutoDiff<N, T> &b)
   {
       T denominator = b.value()*b.value(); //out of the for becuase it is the same for all i
       // so the for does not have to do the same operation each iteration

       AutoDiff<N, T> result(a.value()/b.value());

       // we use the math rul (a/b)'=((a'*b-a*b')/b^2)
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = (a.deriv()[i]*b.value() - a.value()*b.deriv()[i])/ denominator;
       return result;
   }


  // 2. Add some more functions (cos, exp, log, …) for the AutoDiff class.

  // Cos Fucntion 
   template <size_t N, typename T = double>
   AutoDiff<N, T> cos(const AutoDiff<N, T>&a)
   {
       AutoDiff<N, T> result(cos(a.value())); 
       // cos(x)'=-sin(x)*x' (chain rule in case it  is not only an x)
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = -sin(a.value())*a.deriv()[i];
       return result;
   }

  // Exp Fucntion
   using std::exp;

   template <size_t N, typename T = double>
   AutoDiff<N, T> exp(const AutoDiff<N, T>&a)
   {
       T exp_val = exp(a.value()); // out of the loop because it is the same always and this way
       // we don't need to compute it each iteration (to be more efficient)

       AutoDiff<N, T> result(exp_val); 
       // exp(x)'=exp(x)*x'
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = exp_val*a.deriv()[i];
       return result;
   }
   
   // Log function 
   using std::log;

   template <size_t N, typename T = double>
   AutoDiff<N, T> log(const AutoDiff<N, T>&a)
   {
       AutoDiff<N, T> result(log(a.value())); 
       // log(x)'=(1/x)*x' = x'/x
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = a.deriv()[i]/a.value();
       return result;
   }



} // namespace ASC_ode

#endif