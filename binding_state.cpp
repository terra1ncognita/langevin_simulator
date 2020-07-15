#include "binding_state.h"

BindingEvent detect_changes(const bitmask::bitmask<BindingState>& prev, const bitmask::bitmask<BindingState>& curr) {
	return BindingEvent(~prev & curr, prev & ~curr);
}
