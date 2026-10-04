/*
 * newtonUpdaterGeneric.h
 *
 *  Created on: 28.07.2021
 *      Author: Markus Knodel
 */

#ifndef UGCORE_UGBASE_LIB_DISC_OPERATOR_NON_LINEAR_OPERATOR_NEWTON_SOLVER_NEWTONUPDATERGENERIC_H_
#define UGCORE_UGBASE_LIB_DISC_OPERATOR_NON_LINEAR_OPERATOR_NEWTON_SOLVER_NEWTONUPDATERGENERIC_H_

#include <ostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <type_traits>
#include <utility>

#include "common/common.h"

namespace ug
{

template <typename TVector>
class NewtonUpdaterGeneric
{

public:

	virtual ~NewtonUpdaterGeneric() {};


	using vector_type = TVector;

	virtual bool updateSolution( vector_type & sol, vector_type const & corr, bool signNegative = true )
	{
		if( signNegative )
			sol -= corr;
		else
			sol += corr;

		return true;
//		return (*this)(sol, 1, sol, - 1, corr );
	}


	virtual bool updateSolution( vector_type & solNew, number scaleOldSol,  vector_type const & solOld,
					             number scaleCorr, vector_type const & corr )
	{
		VecScaleAdd(solNew, scaleOldSol, solOld, scaleCorr, corr);

		return true;

	}

	virtual bool resetSolution( vector_type & resettedSol, vector_type const & oldSol )
	{
		resettedSol = oldSol;

		// reset not only u, but also internal variables if available!!!!!

		return true;
	}

	virtual bool tellAndFixUpdateEvents( vector_type const & sol )
	{
		return true;
	}
	
	virtual int counter()
	{
		return 0;
	}

};


template <typename TVector>
class NewtonUpdaterProjection : public NewtonUpdaterGeneric<TVector>
{
	public:
		
		// Destructor
		virtual ~NewtonUpdaterProjection() override = default;

		// Constructor
		NewtonUpdaterProjection(): m_count(0), m_numfct(3), u_max(1.1), u_min(1e-05){};
	
		using vector_type = TVector;

		// Override simple updateSolution
		virtual bool updateSolution(vector_type& sol,
									vector_type const& corr,
									bool signNegative = true) override
		{
			// Call base class update
			bool UpdatedSol = NewtonUpdaterGeneric<TVector>::updateSolution(sol, corr, signNegative);

			// Apply projection
			if(UpdatedSol)
				UpdatedSol = UpdatedSol && project(sol);

			return UpdatedSol;
		}

		// Override scaled updateSolution
		virtual bool updateSolution(vector_type& solNew,
									number scaleOldSol,
									vector_type const& solOld,
									number scaleCorr,
									vector_type const& corr) override
		{
			// Call base class update
			bool UpdatedSol = NewtonUpdaterGeneric<TVector>::updateSolution(solNew, scaleOldSol, solOld, scaleCorr, corr);

			// Apply projection
			if(UpdatedSol)
				UpdatedSol = UpdatedSol && project(solNew);

			return UpdatedSol;
		}

		// Override resetSolution
		virtual bool resetSolution(vector_type& resettedSol,
								   vector_type const& oldSol) override
		{
			// Call base class reset
			return NewtonUpdaterGeneric<TVector>::resetSolution(resettedSol, oldSol);
		}
		
		
		// Override tellAndFixUpdateEvents
		virtual bool tellAndFixUpdateEvents( vector_type const & sol ) override
		{
			return NewtonUpdaterGeneric<TVector>::tellAndFixUpdateEvents( sol );

		}
		// Override counter
		virtual int counter( ) override
		{
			return m_count;
		}
	public:
	///	sets variable which will be projected
		void set_projection_fct(int numfct) {m_numfct = numfct;}
	
		void set_projection_mask(const vector_type& mask)
		{
			m_projectionMask = mask;
			m_hasProjectionMask = true;
		}

	///	sets max threshold
		void set_max_threshold(number upper_bound) {u_max = upper_bound;}

	///	sets min threshold
		void set_min_threshold(number lower_bound) {u_min = lower_bound;}
		
	
	protected:
	/// solution copy
		vector_type u_aux;
	/// Main variable index to apply projection
		int m_numfct;
	/// Upper bound
		number u_max;
	/// Lower bound
		number u_min;
		int m_count;
	
		vector_type m_projectionMask;
		bool m_hasProjectionMask = false;
	
	private:

		// Scalar algebra: each entry is a number.
		template <typename T>
		static typename std::enable_if<std::is_arithmetic<T>::value, T&>::type
		projection_value(T& value, int)
		{
			return value;
		}

		// Block algebra: u[k] is a vector/block
		template <typename T>static typename std::enable_if<!std::is_arithmetic<T>::value,decltype(std::declval<T&>()(0,0))>::type
		projection_value(T& value, int fct)
		{
			return value(fct, 0);
		}
	private:
	bool project(vector_type& u)
	{
		const bool scalarAlgebra = std::is_arithmetic<typename vector_type::value_type>::value;

		if(scalarAlgebra)
		{
			if(!m_hasProjectionMask)
				UG_THROW("NewtonUpdaterProjection: CPU1 requires a projection mask.");

			if(m_projectionMask.size() != u.size())
				UG_THROW("NewtonUpdaterProjection: mask and solution sizes differ.");
		}

		if(u_min > u_max)
			UG_THROW("NewtonUpdaterProjection: invalid projection bounds.");

		u_aux.resize(u.size());
		u_aux = u;

		int Count = 0;

		for(size_t k = 0; k < u_aux.size(); ++k)
		{
			// Only scalar algebra uses the mask.
			if(scalarAlgebra)
			{
				if(projection_value(m_projectionMask[k], m_numfct) == 0.0)
					continue;
			}

			// CPU1: scalar entry selected by the mask.
			// Block algebra: component m_numfct inside block k.
			auto& val = projection_value(u_aux[k], m_numfct);

			if(val > u_max || val < u_min)
			{
				val = fmax(u_min, fmin(u_max, val));
				++Count;
			}
		}

		u = u_aux;
		m_count = Count;

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1)
		{
			pcl::ProcessCommunicator com;
			com.allreduce(&Count, &m_count, 1, PCL_DT_INT, PCL_RO_SUM);
		}
	#endif

		return true;
	}
};

}


#endif /* UGCORE_UGBASE_LIB_DISC_OPERATOR_NON_LINEAR_OPERATOR_NEWTON_SOLVER_NEWTONUPDATERGENERIC_H_ */


// muss noch im Newton und line search irgendwie automatisch so initialisiert werden, aber auf User-Wunsch ueberschrieben
// mit Hilfe von gettern und settern, damit entweder der line search oder der Newton selber sich das vom anderen holen koennen
// und dann muss es soweit registriert werden, wie es noetig ist, damit man von aussen zu greifen kann
// und die abgeleiteten Klassen brauchen ggf einen Konstruktor und so weiter
///	Vector type
//typedef typename TAlgebra::vector_type vector_type;

//	using vector_type = TVector; // typename TAlgebra::vector_type;
//	using vec_typ_val_typ = typename vector_type::value_type;

//typedef TVector vector_type; // typename TAlgebra::vector_type;
//typedef typename vector_type::value_type vec_typ_val_typ;
//typedef typename vector_type::value_type value_type;

//	template <typename VecTyp> //, typename ScalTyp >
//	template <typename ScalTyp> //, typename ScalTyp >
//	bool updateNewton( vector_type & solNew, value_type scaleOldSol,  vector_type const & solOld,
//			 	 	   value_type scaleCorr, vector_type const & corr )
//	bool updateNewton( VecTyp & solNew, ScalTyp scaleOldSol,  VecTyp const & solOld,
//			 	 	   ScalTyp scaleCorr, VecTyp const & corr )
//	bool updateNewton( VecTyp & solNew, typename VecTyp::value_type scaleOldSol,  VecTyp const & solOld,
//					   typename VecTyp::value_type scaleCorr, VecTyp const & corr )
//	bool updateNewton( vector_type & solNew, ScalTyp scaleOldSol,  vector_type const & solOld,
//					   ScalTyp scaleCorr, vector_type const & corr )
//		static_assert( std::is_same<ScalTyp, double>::value, "is ScalType double" );
//
//		static_assert( !std::is_same<ScalTyp, number>::value, "! is ScalType number" );
//
//		static_assert( std::is_same<ScalTyp, value_type>::value, "is it value_type" );
//
//		static_assert( std::is_same<double, value_type>::value, "is value_type double" );
//
//		static_assert( std::is_same<number, value_type>::value, "is value_type number" );

		// 	Standard: try on line u := u - lambda*p
		//	VecScaleAdd(u, 1.0, u, (-1)*lambda, p);

//		typedef typename vector_type::size_type size_type;
////	    using size_type = typename vector_type::size_type;
//
//		for(size_type i = 0; i < solNew.size(); ++i)
//		{
//			solNew[i] = scaleOldSol * solOld[i] + scaleCorr * corr[i];
//		}
