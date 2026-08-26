#ifndef WARLOCK_OBJECT_HPP
#define WARLOCK_OBJECT_HPP
#include <cstdint>

namespace framework {
    class Object {
    public:
        Object(double ts = 0, uint64_t ev_id = 0) : timestamp_(ts), event_id_(ev_id) {}
        virtual ~Object() = default;
        
        double timestamp() const { return timestamp_; }
        void setTimestamp(double ts) { timestamp_ = ts; }
        
        uint64_t eventID() const { return event_id_; }
        void setEventID(uint64_t ev_id) { event_id_ = ev_id; }
        
    protected:
        double timestamp_;
        uint64_t event_id_;
    };
}
#endif