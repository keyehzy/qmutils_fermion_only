#pragma once

#include <list>
#include <optional>
#include <unordered_map>
#include <utility>

namespace qmutils {

template <typename Key, typename Value>
class LRUCache {
 public:
  using key_type = Key;
  using value_type = Value;
  using list_type = std::list<std::pair<const Key*, Value>>;
  using map_type = std::unordered_map<Key, typename list_type::iterator>;

  explicit LRUCache(size_t capacity)
      : capacity_(capacity), cache_map_(std::max(capacity, size_t{16})) {
    cache_map_.reserve(capacity);
  }

  void put(const key_type& key, value_type&& value) {
    auto it = cache_map_.find(key);
    if (it != cache_map_.end()) {
      it->second->second = std::move(value);
      cache_list_.splice(cache_list_.begin(), cache_list_, it->second);
      return;
    }

    if (cache_list_.size() >= capacity_) {
      const auto& last_key = *cache_list_.back().first;
      cache_map_.erase(last_key);
      cache_list_.pop_back();
    }

    cache_list_.emplace_front(
        &cache_map_.emplace(key, typename list_type::iterator{}).first->first,
        std::move(value));
    cache_map_[key] = cache_list_.begin();
  }

  void put(const key_type& key, const value_type& value) {
    put(key, value_type(value));
  }

  std::optional<std::reference_wrapper<const value_type>> get(
      const key_type& key) {
    auto it = cache_map_.find(key);
    if (it == cache_map_.end()) {
      return std::nullopt;
    }
    cache_list_.splice(cache_list_.begin(), cache_list_, it->second);
    return std::cref(it->second->second);
  }

  void clear() noexcept {
    cache_map_.clear();
    cache_list_.clear();
  }

  size_t size() const noexcept { return cache_map_.size(); }
  size_t capacity() const noexcept { return capacity_; }

 private:
  const size_t capacity_;
  list_type cache_list_;
  map_type cache_map_;
};

}  // namespace qmutils
