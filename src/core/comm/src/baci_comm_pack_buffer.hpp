/*---------------------------------------------------------------------*/
/*! \file

\brief Data packing for sending over MPI

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_COMM_PACK_BUFFER_HPP
#define BACI_COMM_PACK_BUFFER_HPP

#include "baci_config.hpp"

#include <cstring>
#include <vector>

BACI_NAMESPACE_OPEN

namespace CORE::COMM
{
  class PackBuffer
  {
    friend class SizeMarker;

   public:
    class SizeMarker
    {
     public:
      SizeMarker(PackBuffer& data) : data_(data), oldsize_(0) {}

      ~SizeMarker()
      {
        // set actual object size
        data_.SetObjectSize(oldsize_);
      }

      void Insert()
      {
        // add dummy object size, will be filled later
        int size = 0;
        data_.AddtoPack(size);

        // remember current data size
        oldsize_ = data_().size();
      }

     private:
      PackBuffer& data_;
      std::size_t oldsize_;
    };

    PackBuffer() : size_(0), grow_(true) {}

    void StartPacking()
    {
      grow_ = false;
      buf_.reserve(size_);
    }

    std::vector<char>& operator()() { return buf_; }

    const std::vector<char>& operator()() const { return buf_; }

    /// add POD object
    template <typename kind>
    void AddtoPack(const kind& stuff)
    {
      std::size_t osize = sizeof(kind);
      if (grow_)
      {
        size_ += osize;
      }
      else
      {
        std::size_t oldsize = buf_.size();
        buf_.resize(oldsize + osize);
        std::memcpy(&buf_[oldsize], &stuff, osize);
      }
    }

    /// add array of POD objects
    template <typename kind>
    void AddtoPack(const kind* stuff, std::size_t stuffsize)
    {
      if (grow_)
      {
        size_ += stuffsize;
      }
      else
      {
        std::size_t oldsize = buf_.size();
        buf_.resize(oldsize + stuffsize);
        std::memcpy(&buf_[oldsize], stuff, stuffsize);
      }
    }

   private:
    /// set size of a ParObject after it has been inserted
    void SetObjectSize(std::size_t oldsize)
    {
      if (not grow_)
      {
        int osize = buf_.size() - oldsize;
        std::memcpy(&buf_[oldsize - sizeof(int)], &osize, sizeof(int));
      }
    }

    std::vector<char> buf_;
    std::size_t size_;
    bool grow_;
  };
}  // namespace CORE::COMM

BACI_NAMESPACE_CLOSE

#endif
