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

  using list_type = std::list<std::pair<Key, Value>>;
  using map_type = std::unordered_map<Key, typename list_type::iterator>;

  explicit LRUCache(size_t capacity) : capacity_(capacity) {
    cache_map_.reserve(std::max<size_t>(64, capacity));
  }

  template <typename K, typename V>
  void put(K&& key, V&& value) {
    auto it = cache_map_.find(key);
    if (it != cache_map_.end()) {
      it->second->second = std::forward<V>(value);
      cache_list_.splice(cache_list_.begin(), cache_list_, it->second);
      return;
    }

    if (capacity_ == 0) {
      return;
    }
    if (cache_list_.size() >= capacity_) {
      auto& least_recent_pair = cache_list_.back();
      cache_map_.erase(least_recent_pair.first);
      cache_list_.pop_back();
    }

    cache_list_.emplace_front(std::forward<K>(key), std::forward<V>(value));
    cache_map_[cache_list_.front().first] = cache_list_.begin();
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

  std::optional<std::reference_wrapper<value_type>> get_mutable(
      const key_type& key) {
    auto it = cache_map_.find(key);
    if (it == cache_map_.end()) {
      return std::nullopt;
    }
    cache_list_.splice(cache_list_.begin(), cache_list_, it->second);
    return std::ref(it->second->second);
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
