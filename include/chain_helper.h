#pragma once

#include <TBranch.h>
#include <TChain.h>
#include <concepts>
#include <type_traits>
#include <memory>
#include <string>
#include <vector>

template <typename... branch_type> class root_chain {
private:
    std::unique_ptr<TChain> chain{};
    std::tuple<branch_type *...> data;
    std::array<TBranch *, sizeof...(branch_type)> b_add;

    template <class T, T I, T... J> inline void do_delete(std::integer_sequence<T, I, J...>) {
        do_delete(std::integer_sequence<std::size_t, I>{});
        do_delete(std::integer_sequence<std::size_t, J...>{});
    }

    template <class T, T I> inline void do_delete(std::integer_sequence<T, I>) {
        typedef typename std::tuple_element<I, std::tuple<branch_type...>>::type thistype;
        if (std::is_trivially_copyable<thistype>::value)
            std::default_delete<thistype>{}(std::get<I>(data));
        // simliarI to std::is_array<thistype>::value?delete [] std::get<>(data):delete std::get<I>(data);
    }

    template <class T, T I, T... J> inline void do_setaddress(std::integer_sequence<T, I, J...>, const std::array<const char *, sizeof...(branch_type)> names) {
        do_setaddress(std::integer_sequence<std::size_t, I>{}, names);
        do_setaddress(std::integer_sequence<std::size_t, J...>{}, names);
    }

    template <class T, T I> inline void do_setaddress(std::integer_sequence<T, I>, const std::array<const char *, sizeof...(branch_type)> names) {
        b_add[I] = chain->GetBranch(std::get<I>(names));
        assert(b_add[I]);
        typedef typename std::tuple_element<I, std::tuple<branch_type...>>::type thistype;
        if (std::is_trivially_copyable<thistype>::value) {
            // normal case, assign space and SetBranchAddress, the ugly expression
            // is for properly dealing with arrays, the assigned space is then freed
            // by me (not ROOT)
            std::get<I>(data) = reinterpret_cast<thistype *>(new thistype);
            chain->SetBranchAddress(std::get<I>(names), std::get<I>(data), &b_add[I]);
            // following method is broken, don't know why:
            // b_add[I]->SetAddress(std::get<I>(data));
        } else {
            // ok, seems for non-trivially copyable containers, ROOT will only accept
            // T** pointer, and assign space on its own, free on its own
            std::get<I>(data) = nullptr;
            chain->SetBranchAddress(std::get<I>(names), &std::get<I>(data), &b_add[I]);
        }
    }

    template <typename T, T... IDs> std::tuple<const branch_type &...> get_elements(std::integer_sequence<T, IDs...>, std::size_t id) {
        auto curr = chain->LoadTree(id);
        if (curr < 0) {
            throw;
        }
        for (auto &i : b_add)
            i->GetEntry(curr);
        // chain->GetEntry(curr);
        return {*std::get<IDs>(data)...};
    }

public:
    template <typename T> 
    requires std::is_same_v<std::remove_cvref_t<typename std::remove_cvref_t<T>::value_type>, std::string>
    root_chain(T &&file_list, const char *tree_name, const std::array<const char *, sizeof...(branch_type)> names) {
        chain = std::make_unique<TChain>(tree_name);
        for (const auto &i : file_list) {
            chain->Add(i.c_str());
        }
        do_setaddress(std::make_integer_sequence<std::size_t, sizeof...(branch_type)>{}, names);
    }

    root_chain(const root_chain &) = delete;

    ~root_chain() {
        chain->ResetBranchAddresses();
        do_delete(std::make_integer_sequence<std::size_t, sizeof...(branch_type)>{});
        // cleaning of chain is done by unique_ptr
        // the chain will then clean up TBranches
    }

    decltype(auto) get_elements(std::size_t id) { return get_elements(std::make_integer_sequence<std::size_t, sizeof...(branch_type)>{}, id); }

    decltype(auto) get_entries() { return chain->GetEntries(); }

    // used by range-based-for
    class iter {
    private:
        std::size_t num;
        root_chain &se;

    public:
        iter(int m_num, root_chain &root_chain_self) : num(m_num), se(root_chain_self) {}
        iter &operator++() {
            num++;
            return *this;
        }
        bool operator!=(const iter &other) const { return num != other.num; }
        decltype(auto) operator*() { return se.get_elements(num); }
    };

    decltype(auto) operator[](std::size_t id) { return get_elements(id); }

    decltype(auto) begin() { return iter(0, *this); }
    decltype(auto) end() { return iter(get_entries(), *this); }
};
