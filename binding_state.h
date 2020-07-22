#pragma once
#include "bitmask.hpp"

enum class BindingState {
	Free = 0,
	FirstSiteBound = 1 << 0,
	SecondSiteBound = 1 << 1,
	NewWell = 1 << 2,

	_bitmask_max_element = NewWell,
};

BITMASK_DEFINE(BindingState);

//
//template<typename T, typename = std::enable_if<std::is_enum<T>::value, T>::type>
//class State
//{
//	T val_;
//public:
//	constexpr State(T val) : val_(val) {}
//	constexpr operator T() const { return val_; }
//	explicit operator bool() const {
//		return static_cast<std::underlying_type<T>::type>(val_) != 0;
//	}
//	constexpr bool contains(const T other) {
//		return bool(val_ & other);
//	}
//};
//
//// LOGICAL AND
//template<typename T, typename = std::enable_if<std::is_enum<T>::value, T>::type>
//constexpr State<T> operator &(const T x, const T y) {
//	return State<T>(static_cast<typename std::underlying_type<T>::type>(x) &
//		static_cast<typename std::underlying_type<T>::type>(y));
//}
//template<typename T, typename = std::enable_if<std::is_enum<T>::value, T>::type>
//constexpr State<T> operator &(const State<T> x, const State<T> y) {
//	return State<T>(static_cast<typename std::underlying_type<T>::type>(x) &
//		static_cast<typename std::underlying_type<T>::type>(y));
//}
//
//
//// LOGICAL OR
//template<typename T, typename = std::enable_if<std::is_enum<T>::value, T>::type>
//constexpr State<T> operator |(const T x, const T y) {
//	return State<T>(static_cast<typename std::underlying_type<T>::type>(x) |
//		static_cast<typename std::underlying_type<T>::type>(y));
//}
//
//// LOGICAL NOT
//template<typename T, typename = std::enable_if<std::is_enum<T>::value, T>::type>
//constexpr State<T> operator ~(const T x) {
//	return State<T>(~static_cast<typename std::underlying_type<T>::type>(x));
//}
//template<typename T, typename = std::enable_if<std::is_enum<T>::value, T>::type>
//constexpr State<T> operator ~(const State<T> x) {
//	return State<T>(~(x()));
//}
//
//// OUTPUT
//template<typename T, typename = std::enable_if<std::is_enum<T>::value, T>::type>
//std::ostream& operator <<(std::ostream& out, const T x) {
//	out << (static_cast<typename std::underlying_type<T>::type>(x));
//	return out;
//}
//

class BindingEvent {
public:
	bitmask::bitmask<BindingState> binding, release;
	BindingEvent(const bitmask::bitmask<BindingState>& _binding = 0, const bitmask::bitmask<BindingState>& _release = 0) : binding(_binding), release(_release) {}

	bool any_change() {
		return bool(binding | release);
	}
};

template<typename Char, typename Traits>
auto operator<<(std::basic_ostream<Char, Traits> &os, const BindingEvent& x)
-> std::basic_ostream<Char, Traits>& {
	os << x.binding.bits() << x.release.bits();
	return os;
}

BindingEvent detect_changes(const bitmask::bitmask<BindingState>& prev, const bitmask::bitmask<BindingState>& curr);
