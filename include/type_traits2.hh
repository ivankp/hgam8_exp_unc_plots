#ifndef IVANP_TYPE_TRAITS
#define IVANP_TYPE_TRAITS

template <typename,typename> struct change_value_type { };
template <typename Tnew, typename Told, typename Alloc>
struct change_value_type<std::vector<Told,Alloc>,Tnew> {
  using type = std::vector<Tnew>;
};
template <typename Tnew, typename Told, size_t N>
struct change_value_type<std::array<Told,N>,Tnew> {
  using type = std::array<Tnew,N>;
};

template <typename T>
struct transpose_container_type {
private:
  using outer = T;
  using inner = typename outer::value_type;
  using value_type = typename inner::value_type;
public:
  using type = typename change_value_type<inner,
    typename change_value_type<outer,value_type>::type >::type;
};
template <typename T>
using transpose_container_t = typename transpose_container_type<T>::type;

#endif
