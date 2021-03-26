#pragma once
#include <immintrin.h>
#include <random>
#include <type_traits>

namespace PackedCSparse {
	template <typename T, size_t k>
	struct BaseImpl
	{
		T x[k];

		BaseImpl() {};

		BaseImpl(const T& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] = rhs;
		}

		BaseImpl operator+(const BaseImpl& rhs) const
		{
			BaseImpl lhs;
			for (size_t i = 0; i < k; i++)
				lhs.x[i] = x[i] + rhs.x[i];
			return lhs;
		}

		BaseImpl operator-(const BaseImpl& rhs) const
		{
			BaseImpl lhs;
			for (size_t i = 0; i < k; i++)
				lhs.x[i] = x[i] - rhs.x[i];
			return lhs;
		}

		BaseImpl operator*(const BaseImpl& rhs) const
		{
			BaseImpl lhs;
			for (size_t i = 0; i < k; i++)
				lhs.x[i] = x[i] * rhs.x[i];
			return lhs;
		}

		BaseImpl operator/(const BaseImpl& rhs) const
		{
			BaseImpl lhs;
			for (size_t i = 0; i < k; i++)
				lhs.x[i] = x[i] / rhs.x[i];
			return lhs;
		}

		BaseImpl& operator+=(const BaseImpl& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] += rhs.x[i];
			return *this;
		}

		BaseImpl& operator-=(const BaseImpl& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] -= rhs.x[i];
			return *this;
		}

		BaseImpl& operator*=(const BaseImpl& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] *= rhs.x[i];
			return *this;
		}

		BaseImpl& operator/=(const BaseImpl& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] /= rhs.x[i];
			return *this;
		}

		explicit operator bool() const
		{
			bool ret = false;
			for (size_t i = 0; i < k; i++)
				ret = ret || bool(x[i]);
			return ret;
		}

		static T get(const BaseImpl& a, size_t index)
		{
			return a.x[index];
		}

		static void set(BaseImpl& a, size_t index, const T& value)
		{
			a.x[index] = value;
		}

		static BaseImpl abs(BaseImpl& x)
		{
			BaseImpl out;
			for (size_t i = 0; i < k; i++)
				out.x[i] = abs(x[i]);
			return out;
		}

		static BaseImpl log(BaseImpl& x)
		{
			BaseImpl out;
			for (size_t i = 0; i < k; i++)
				out.x[i] = log(x[i]);
			return out;
		}

		static void fmadd(BaseImpl& a, const BaseImpl& b, const BaseImpl& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] += b.x[i] * c.x[i];
		}

		static void fnmadd(BaseImpl& a, const BaseImpl& b, const BaseImpl& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] -= b.x[i] * c.x[i];
		}

		static void fmadd(BaseImpl& a, const BaseImpl& b, const T& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] += b.x[i] * c;
		}

		static void fnmadd(BaseImpl& a, const BaseImpl& b, const T& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] -= b.x[i] * c;
		}

		static BaseImpl clipped_sqrt(const BaseImpl& a, const T nonpos_output)
		{
			BaseImpl out;
			for (size_t i = 0; i < k; i++)
			{
				T r = a.x[i];
				if (r > 0)
					out.x[i] = sqrt(r);
				else
					out.x[i] = nonpos_output;
			}
			return out;
		}

		static BaseImpl sign(std::mt19937_64& gen)
		{
			BaseImpl out;
			unsigned long long seed = gen();
			for (size_t i = 0; i < k; i++)
			{
				out.x[i] = T((2 * ((seed >> i) & 1)) - 1.0);
				if ((i & 63) == 63) seed = gen();
			}
			return out;
		}

	};

	template <typename T>
	struct BaseScalarImpl
	{
		static T get(const T& x, size_t index)
		{
			return x;
		}

		static void set(T& x, size_t index, T& value)
		{
			x = value;
		}

		static T abs(const T &x)
		{
			return ::abs(x);
		}

		static T log(const T &x)
		{
			return ::log(x);
		}

		static void fmadd(T& a, const T& b, const T& c)
		{
			a += b * c;
		}

		static void fnmadd(T& a, const T& b, const T& c)
		{
			a -= b * c;
		}

		static T clipped_sqrt(const T& x, const T& nonpos_output)
		{
			if (x > 0)
				return sqrt(x);
			else
				return nonpos_output;
		}

		static T sign(std::mt19937_64& gen)
		{
			unsigned long long seed = gen();
			return T((2 * (seed & 1)) - 1.0);
		}
	};

	template<size_t k>
	struct m256dArray
	{
		__m256d x[k];

		m256dArray() {};

		m256dArray(const double rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] = _mm256_set1_pd(rhs);
		}

		template<size_t k2>
		m256dArray(const m256dArray<k2>& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] = rhs.x[i % k2];
		}

		m256dArray operator+(const m256dArray& rhs) const
		{
			m256dArray out;
			for (size_t i = 0; i < k; i++)
				out.x[i] = _mm256_add_pd(x[i], rhs.x[i]);
			return out;
		}

		m256dArray operator-(const m256dArray& rhs) const
		{
			m256dArray out;
			for (size_t i = 0; i < k; i++)
				out.x[i] = _mm256_sub_pd(x[i], rhs.x[i]);
			return out;
		}

		m256dArray operator*(const m256dArray& rhs) const
		{
			m256dArray out;
			for (size_t i = 0; i < k; i++)
				out.x[i] = _mm256_mul_pd(x[i], rhs.x[i]);
			return out;
		}

		m256dArray operator/(const m256dArray& rhs) const
		{
			m256dArray out;
			for (size_t i = 0; i < k; i++)
				out.x[i] = _mm256_div_pd(x[i], rhs.x[i]);
			return out;
		}

		m256dArray& operator+=(const m256dArray& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] = _mm256_add_pd(x[i], rhs.x[i]);
			return *this;
		}

		m256dArray& operator-=(const m256dArray& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] = _mm256_sub_pd(x[i], rhs.x[i]);
			return *this;
		}

		m256dArray& operator*=(const m256dArray& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] = _mm256_mul_pd(x[i], rhs.x[i]);
			return *this;
		}

		m256dArray& operator/=(const m256dArray& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] = _mm256_div_pd(x[i], rhs.x[i]);
			return *this;
		}

		explicit operator bool() const
		{
			bool ret = false;
			__m256d z = _mm256_set1_pd(0.0);
			for (size_t i = 0; i < k; i++)
			{
				__m256d c = _mm256_cmp_pd(x[i], z, _CMP_EQ_OQ);
				ret = ret || (_mm256_movemask_pd(c) != 0xf);
			}
			return ret;
		}

		static double get(const m256dArray& x, size_t index)
		{
			double y[4];
			_mm256_store_pd(y, x.x[index / 4]);
			return y[index & 3];
		}

		static void set(m256dArray& x, size_t index, double value)
		{
			__m256d v = _mm256_broadcast_sd(&value);
			switch (index & 3)
			{
			case 0:  x.x[index / 4] = _mm256_blend_pd(x.x[index / 4], v, 1); break;
			case 1:  x.x[index / 4] = _mm256_blend_pd(x.x[index / 4], v, 2); break;
			case 2:  x.x[index / 4] = _mm256_blend_pd(x.x[index / 4], v, 4); break;
			default: x.x[index / 4] = _mm256_blend_pd(x.x[index / 4], v, 8); break;
			}
		}

		static m256dArray abs(const m256dArray& x)
		{
			const __m256d mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));

			m256dArray out;
			for (size_t i = 0; i < k; i++)
				out.x[i] = _mm256_and_pd(x.x[i], mask);
			return out;
		}

		static m256dArray log(const m256dArray& x)
		{
            // gcc does not support _mm256_log_pd
            // Do it sequentially instead
            
			//m256dArray out;
			//for (size_t i = 0; i < k; i++)
			//	out.x[i] = _mm256_log_pd(x.x[i]);
            
            m256dArray out;
            for (size_t i = 0; i < 4*k; i++)
                set(out, i, std::log(get(x,i)));
			return out;
		}

		static void fmadd(m256dArray& a, const m256dArray& b, const double& c)
		{
			auto cx = _mm256_set1_pd(c);
			for (size_t i = 0; i < k; i++)
				a.x[i] = _mm256_fmadd_pd(b.x[i], cx, a.x[i]);
		}

		static void fnmadd(m256dArray& a, const m256dArray& b, const double& c)
		{
			auto cx = _mm256_set1_pd(c);
			for (size_t i = 0; i < k; i++)
				a.x[i] = _mm256_fnmadd_pd(b.x[i], cx, a.x[i]);
		}

		static void fmadd(m256dArray& a, const m256dArray& b, const m256dArray& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] = _mm256_fmadd_pd(b.x[i], c.x[i], a.x[i]);
		}

		static void fnmadd(m256dArray& a, const m256dArray& b, const m256dArray& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] = _mm256_fnmadd_pd(b.x[i], c.x[i], a.x[i]);
		}
		
		static m256dArray clipped_sqrt(const m256dArray& x, const double nonpos_output)
		{
			m256dArray out;

			const __m256d large = { nonpos_output, nonpos_output, nonpos_output, nonpos_output };
			const __m256d zero = _mm256_setzero_pd();
			for (size_t i = 0; i < k; i++)
			{
				__m256d xi = x.x[i];
				__m256d mask = _mm256_cmp_pd(xi, zero, _CMP_LE_OS); // mask = (rhs.x[i]<= 0) ? -1 : 0
				out.x[i] = _mm256_blendv_pd(_mm256_sqrt_pd(xi), large, mask);
			}
			return out;
		}

		static m256dArray sign(std::mt19937_64& gen)
		{
			m256dArray out;
			const __m256i bits = _mm256_set_epi64x(1, 2, 4, 8);
			const __m256d zero = _mm256_setzero_pd();
			const __m256d pos = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
			const __m256d neg = _mm256_set_pd(-1.0, -1.0, -1.0, -1.0);

			unsigned long long seed = gen();
			for (size_t i = 0; i < k; i++)
			{
				__m256i s = _mm256_set1_epi64x((seed >> (4 * i)) & 15);
				__m256i xi = _mm256_and_si256(s, bits);
				__m256d x = _mm256_castsi256_pd(xi);
				__m256d mask = _mm256_cmp_pd(x, zero, _CMP_EQ_OQ); // mask = (rhs.x[i] == 0) ? -1 : 0
				out.x[i] = _mm256_blendv_pd(pos, neg, mask);
				if ((i & 63) == 63) seed = gen();
			}
			return out;
		}
	};


	template <typename T, size_t k>
	struct FloatTypeSelector
	{
		using type = typename std::conditional<k == 1, T, BaseImpl<T, k>>::type;
		using funcImpl = typename std::conditional<k == 1, BaseScalarImpl<T>, BaseImpl<T, k>>::type;
	};

	template <size_t k>
	struct FloatTypeSelector<double, k>
	{
		static_assert(k == 1 || k % 4 == 0, "Array<double,k> assumes k = 1 or a multiple of 4");
		using type = typename std::conditional< k == 1, double, m256dArray<k / 4>>::type;
		using funcImpl = typename std::conditional< k == 1, BaseScalarImpl<double>, m256dArray<k / 4>>::type;
	};

	template <size_t k, size_t l>
	struct FloatTypeSelector<m256dArray<k>, l>
	{
		using type = m256dArray<k* l>;
		using funcImpl = m256dArray<k* l>;
	};

	template <typename T, size_t k, size_t l>
	struct FloatTypeSelector<BaseImpl<T, k>, l>
	{
		using type = BaseImpl<T, k* l>;
		using funcImpl = BaseImpl<T, k* l>;
	};

	template<typename T, size_t k = 1>
	using FloatArray = typename FloatTypeSelector<T, k>::type;

	template<typename T, size_t k = 1>
	using FloatArrayFunc = typename FloatTypeSelector<T, k>::funcImpl;

	template<typename T>
	auto get(const T& a, size_t index) -> decltype(FloatArrayFunc<T>::get(a, index))
	{
		return FloatArrayFunc<T>::get(a, index);
	}

	template<typename T1, typename T2>
	void set(T1& a, size_t index, T2 value)
	{
		FloatArrayFunc<T1>::set(a, index, value);
	}

	template<typename T1, typename T2, typename T3>
	void fmadd(T1& a, const T2& b, const T3& c)
	{
		FloatArrayFunc<T1>::fmadd(a, b, c);
	}

	template<typename T1, typename T2, typename T3>
	void fnmadd(T1& a, const T2& b, const T3& c)
	{
		FloatArrayFunc<T1>::fnmadd(a, b, c);
	}

	template<typename T1, typename T2>
	T1 clipped_sqrt(const T1& a, const T2 b)
	{
		return FloatArrayFunc<T1>::clipped_sqrt(a, b);
	}

	template<typename T>
	T abs(const T& a)
	{
		return FloatArrayFunc<T>::abs(a);
	}

	template<typename T>
	T log(const T& a)
	{
		return FloatArrayFunc<T>::log(a);
	}

	template<typename T>
	T sign(std::mt19937_64& gen)
	{
		return FloatArrayFunc<T>::sign(gen);
	}
}