
#if 1
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <array>
#include <deque>
#include <algorithm>
#include <utility>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <numeric>
#include <assert.h>
#include <bitset>
#include <list>
#include <cmath>
#include <execution>
#include <random>
#include "boost/range/adaptor/indexed.hpp"

#define all_range(C) std::begin(C), std::end(C)

constexpr int LENGTH_SUM_MAX = 10000;
constexpr int LENGTH_MAX = 500;

template<typename T = int_fast32_t>
struct Range {
private:
	using This_T = Range<T>;
	using value_type = T;
	value_type b, e;
public:
	constexpr Range(value_type end)noexcept :b(0), e(end) {}
	constexpr Range(value_type begin, value_type end)noexcept :b(begin), e(end) {}
	struct iterator {
	public:
		using value_type = This_T::value_type;
		using pointer = value_type*;
		using difference_type = int;
		using reference = value_type&;
		using const_reference = const value_type&;
		using iterator_category = std::forward_iterator_tag;
	private:
		value_type i;
		friend This_T;
		constexpr iterator(value_type init)noexcept :i(init) {}
	public:
		constexpr iterator()noexcept :i() {}

		constexpr iterator& operator++()noexcept { ++i; return *this; }
		constexpr iterator operator++(int)noexcept { auto tmp = i; ++i; return iterator(tmp); }
		constexpr const_reference operator*() const noexcept { return i; }

		constexpr bool operator==(const iterator& b) const noexcept { return **this == *b; }
		constexpr bool operator!=(const iterator& b) const noexcept { return !((*this) == b); }
	};

	constexpr iterator begin()const noexcept { return iterator(b); }
	constexpr iterator end()const noexcept { return iterator(e); }
};

template<typename Iter>
class indexed_iterator {
public:
	using value_type = std::pair<int, typename Iter::value_type>;
	using pointer = std::pair<int, typename Iter::pointer>;
	using difference_type = int;
	using reference = std::pair<int, typename Iter::reference>;
	using const_reference = std::pair<int, const typename Iter::value_type&>;
	using iterator_category = std::forward_iterator_tag;
private:
	int i = 0;
	Iter iter;
public:
	indexed_iterator()noexcept :i(), iter() {}
	indexed_iterator(Iter iter)noexcept :i(), iter(std::move(iter)) {}

	indexed_iterator& operator++()noexcept { ++i; ++iter; return *this; }
	indexed_iterator operator++(int)noexcept { auto tmp = *this; ++* this; return iterator(tmp); }
	const_reference operator*() const noexcept { return const_reference(i, *iter); }

	bool operator==(const indexed_iterator& b) const noexcept { return iter == b.iter; }
	bool operator!=(const indexed_iterator& b) const noexcept { return !((*this) == b); }
};
template<typename Iter>
indexed_iterator<std::unwrap_reference_t<Iter>> make_indexed_iterator(Iter&& iter) {
	return { std::forward<Iter>(iter) };
}

enum class Type :int { A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, B, Z, X, GAP, NUM };
struct Profile {
	std::array<int32_t, static_cast<int>(Type::NUM)> typenum;//それぞれのtypeの個数
};
inline constexpr int_fast32_t calc_score(const Type& seq1, const Type& seq2) {

	constexpr int_fast32_t table[static_cast<int>(Type::NUM)][static_cast<int>(Type::NUM)] = {
		//  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   - 
		  { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4}, // A 
		  {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4}, // R 
		  {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4}, // N 
		  {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4}, // D 
		  { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4}, // C 
		  {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4}, // Q 
		  {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4}, // E 
		  { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4}, // G 
		  {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4}, // H 
		  {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4}, // I 
		  {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4}, // L 
		  {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4}, // K 
		  {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4}, // M 
		  {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4}, // F 
		  {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4}, // P 
		  { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4}, // S 
		  { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4}, // T 
		  {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4}, // W 
		  {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4}, // Y 
		  { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4}, // V 
		  {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4}, // B 
		  {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4}, // Z 
		  { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4}, // X 
		  {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}, // - 
	};
	return table[static_cast<int>(seq1)][static_cast<int>(seq2)];
}
constexpr int_fast32_t calc_score(const Type& seq, const Profile& prof) {
	int_fast32_t score = 0;
	for (size_t i = 0; i < static_cast<size_t>(Type::NUM); i++) {
		score += prof.typenum[i] * calc_score(seq, static_cast<Type>(i));
	}
	return score;
}
constexpr int_fast32_t calc_profile_score(const Profile& prof) {
	int_fast32_t score = 0;
	for (size_t i = 0; i < static_cast<size_t>(Type::NUM); i++) {
		if (prof.typenum[i] <= 0) { continue; }
		score += prof.typenum[i] * (prof.typenum[i] - 1) * calc_score(static_cast<Type>(i), static_cast<Type>(i)) / 2;
		for (size_t j = i + 1; j < static_cast<size_t>(Type::NUM); j++) {
			score += prof.typenum[i] * prof.typenum[j] * calc_score(static_cast<Type>(i), static_cast<Type>(j));
		}
	}
	return score;
}
int_fast32_t calc_profiles_score(const std::vector<Profile>& profiles) {
	int_fast32_t score = 0;
	for (auto& prof : profiles) {
		score += calc_profile_score(prof);
	}
	return score;
}

struct AlignData {
	//メモリ節約の為、特殊な保持をする
	int_fast32_t score = 0;
	//std::shared_ptr<std::map<size_t, int32_t>> globalGAP;//全てのものについて[n]番目にglobalGAP[n]個ギャップを入れる
	//std::vector<std::shared_ptr<std::map<size_t, int32_t>>> globalGAP_diff;//globalGAPを無効化する分
	//std::vector<std::shared_ptr<std::vector<int32_t>>> localGAP;//[i]番目の文字は、文字[j]の直前にaddLocalGAP[i][j]個のGAPを入れる
	std::vector<std::vector<bool>> localGAP;//1ならGAP
};
inline bool operator<(const AlignData& l, const AlignData& r) {
	return l.score < r.score;
}

std::vector<std::vector<Type>> targetSequences;

std::vector<Type> align_to_array(const AlignData& align, size_t sequence_index)
{
	const auto& si = sequence_index;
	const auto& target = targetSequences[si];
	//auto& globalGAP = *align.globalGAP;
	//const auto& globalGAP_diff = *align.globalGAP_diff[si];
	const auto& localGAP = align.localGAP[si];
	const auto get = [](const std::map<size_t, int32_t>& map, size_t key) {
		auto iter = map.find(key);
		if (iter != map.end()) {
			return iter->second;
		}
		return 0;
	};

	std::vector<Type> res;
	//res.reserve(target.size() * 4 / 3);
	res.reserve(localGAP.size() - 1);
	size_t i = 0;
	//int gaprem = localGAP[0];
	while (true) {

		const size_t n = res.size();
		if (localGAP[n]) {
			res.push_back(Type::GAP);
			continue;
		}
		if (i == target.size()) {
			break;
		}
		res.push_back(target[i]);
		++i;


		////globalGAPの処理
		//auto globalG = get(globalGAP, n) - get(globalGAP_diff, n);
		//if (globalG != 0) {
		//	assert(globalG > 0);
		//	gaprem += globalG;
		//}
		//if (gaprem > 0) {
		//	--gaprem;
		//	res.push_back(Type::GAP);
		//	continue;
		//}

		////文字を進める
		//if (i == target.size()) {
		//	break;
		//}
		//res.push_back(target[i]);
		//++i;

		////localGAPの処理
		//gaprem += localGAP[i];
	}

	return std::move(res);
}
//O(|S||FS|)
void SequenceFixedSequenceAlign(
	int align_sequence, int fixed_sequence, AlignData& output)
{
	const std::vector<Type>& sequence = targetSequences[align_sequence];
	const std::vector<Type>& f_sequence = targetSequences[fixed_sequence];
	static int_fast32_t DP_MEMORY[LENGTH_MAX + 1][LENGTH_MAX + 1];

	size_t si = 0;//sequence_index
	size_t fi = 0;//fixed_sequence_index
	DP_MEMORY[si][fi] = 0;
	for (fi = 1; fi <= f_sequence.size(); ++fi)
	{
		DP_MEMORY[si][fi] = DP_MEMORY[si][fi - 1] + calc_score(Type::GAP, f_sequence[fi - 1]);
	}
	for (si = 1; si <= sequence.size(); ++si)
	{
		//auto& dp_sprev = DP_MEMORY[(si & 1) ^ 1];
		//auto& dp = DP_MEMORY[si & 1];
		auto& dp_sprev = DP_MEMORY[si - 1];
		auto& dp = DP_MEMORY[si];
		fi = si;
		dp[fi] = dp_sprev[fi - 1] + calc_score(sequence[si - 1], f_sequence[fi - 1]);//比較
		for (fi = si + 1; fi <= static_cast<int>(f_sequence.size()); ++fi)
		{
			dp[fi] = std::max({
				dp[fi - 1] + calc_score(Type::GAP, f_sequence[fi - 1]),//sGAP
				dp_sprev[fi - 1] + calc_score(sequence[si - 1], f_sequence[fi - 1]),//比較
				});
		}
	}

	//経路情報
	int32_t s_aligni = sequence.size();//sequence_index
	int32_t fi_aligni = f_sequence.size();//sequence_index
	while (s_aligni > 0 || fi_aligni > 0) {
		auto& dp_sprev = DP_MEMORY[s_aligni - 1 >= 0 ? s_aligni - 1 : 0];
		auto& dp = DP_MEMORY[s_aligni];
		if (s_aligni > 0 && fi_aligni > 0 && dp[fi_aligni] == dp_sprev[fi_aligni - 1] + calc_score(sequence[s_aligni - 1], f_sequence[fi_aligni - 1])) {
			//比較
			--fi_aligni;
			--s_aligni;
			output.localGAP[align_sequence][fi_aligni] = false;
			continue;
		}
		if (fi_aligni > 0 && dp[fi_aligni] == dp[fi_aligni - 1] + calc_score(Type::GAP, f_sequence[fi_aligni - 1])) {
			//sGAP
			--fi_aligni;
			output.localGAP[align_sequence][fi_aligni] = true;
			continue;
		}
		assert(false);
	}

	//DEBUG
#ifdef _DEBUG
	{
		auto r_prof = align_to_array(output, align_sequence);
		if (r_prof.size() != f_sequence.size()) {
			assert(r_prof.size() == f_sequence.size());
		}
	}
#endif
	return;
}
//O(|S||P|)
std::pair<AlignData, std::vector<Profile>> SequenceProfileAlign(
	int align_sequence, const std::vector<Profile>& profile, const AlignData& base_align,
	int_fast32_t(&DP_MEMORY)[LENGTH_MAX + 1][2 * LENGTH_MAX + 1])
{
	const std::vector<Type>& sequence = targetSequences[align_sequence];
	assert(!sequence.empty());
	assert(!profile.empty());
	Profile GapProfile;
	fill(all_range(GapProfile.typenum), 0);
	GapProfile.typenum[static_cast<int>(Type::GAP)] = std::accumulate(all_range(profile[0].typenum), (decltype(profile[0].typenum[0]))0);

	//DP
	int si = 0;//sequence_index
	int pi = 0;//profile_index
	//DP_MEMORY[pi] := sequenceをsi個、profileをpi個使ったときのスコア
	//最終的にはsi==sequence.size()でDP_MEMORY[profile.size()]を取る

	DP_MEMORY[si][pi] = 0;
	si = 0;
	for (pi = 1; pi <= static_cast<int>(profile.size()); ++pi)
	{
		DP_MEMORY[si][pi] = DP_MEMORY[si][pi - 1] + calc_score(Type::GAP, profile[pi - 1]);
	}
	pi = 0;
	for (si = 1; si <= static_cast<int>(sequence.size()); ++si)
	{
		DP_MEMORY[si][pi] = DP_MEMORY[si - 1][pi] + calc_score(sequence[si - 1], GapProfile);
	}
	for (int depth = 2; depth <= static_cast<int>(sequence.size() + profile.size()); depth++)
	{
		Range range(std::max(depth - static_cast<int>(sequence.size()), 1), 1 + std::min(static_cast<int>(profile.size()), depth - 1));

		std::for_each(std::execution::par_unseq, all_range(range), [depth, &DP_MEMORY, &sequence, &profile, &GapProfile](int index) {
			int si = depth - index;
			int pi = index;
			DP_MEMORY[si][pi] = std::max({
				DP_MEMORY[si][pi - 1] + calc_score(Type::GAP, profile[pi - 1]),//sGAP
				//DP_MEMORY[si - 1][pi] + calc_score(sequence[si - 1], GapProfile),//pGAP
				DP_MEMORY[si - 1][pi] + GapProfile.typenum[static_cast<int>(Type::GAP)] * calc_score(sequence[si - 1], Type::GAP),//pGAP
				DP_MEMORY[si - 1][pi - 1] + calc_score(sequence[si - 1], profile[pi - 1]),//比較
				});
		});
	}
	//for (si = 1; si <= static_cast<int>(sequence.size()); ++si)
	//{
	//	//auto& dp_sprev = DP_MEMORY[(si & 1) ^ 1];
	//	//auto& dp = DP_MEMORY[si & 1];
	//	auto& dp_sprev = DP_MEMORY[si - 1];
	//	auto& dp = DP_MEMORY[si];
	//	pi = 0;
	//	dp[pi] = dp_sprev[pi] + calc_score(sequence[si - 1], GapProfile); //pGAP
	//	for (pi = 1; pi <= static_cast<int>(profile.size()); ++pi)
	//	{
	//		dp[pi] = std::max({
	//			dp[pi - 1] + calc_score(Type::GAP, profile[pi - 1]),//sGAP
	//			dp_sprev[pi] + calc_score(sequence[si - 1], GapProfile),//pGAP
	//			dp_sprev[pi - 1] + calc_score(sequence[si - 1], profile[pi - 1]),//比較
	//			});
	//	}
	//}

	//経路情報
	int32_t s_aligni = sequence.size();//sequence_index
	int32_t p_aligni = profile.size();//profile_index
	std::vector<std::pair<int32_t, int32_t>> match;
	match.push_back({ s_aligni , p_aligni });
	while (s_aligni > 0 || p_aligni > 0) {
		auto& dp_sprev = DP_MEMORY[s_aligni - 1 >= 0 ? s_aligni - 1 : 0];
		auto& dp = DP_MEMORY[s_aligni];
		if (p_aligni > 0 && dp[p_aligni] == dp[p_aligni - 1] + calc_score(Type::GAP, profile[p_aligni - 1])) {
			//sGAP
			--p_aligni;
			if (profile[p_aligni].typenum == GapProfile.typenum) {
				continue;
			}
			match.push_back({ s_aligni , p_aligni });
			continue;
		}
		assert(s_aligni > 0);
		if (dp[p_aligni] == dp_sprev[p_aligni] + calc_score(sequence[s_aligni - 1], GapProfile)) {
			//pGAP
			--s_aligni;
			match.push_back({ s_aligni , p_aligni });
			continue;
		}
		assert(p_aligni > 0);
		assert(dp[p_aligni] == dp_sprev[p_aligni - 1] + calc_score(sequence[s_aligni - 1], profile[p_aligni - 1]));
		{
			//比較
			--p_aligni;
			--s_aligni;
			match.push_back({ s_aligni , p_aligni });
			continue;
		}
		//assert(false);
	}
	match.pop_back();//0,0
	std::reverse(all_range(match));


	std::vector<Profile> result_profile;
	result_profile.reserve(match.size());
	s_aligni = match[0].first;//sequence_index
	p_aligni = match[0].second;//profile_index
	AlignData align;
	//align.globalGAP = std::make_unique<decltype(align.globalGAP)::element_type>(*base_align.globalGAP);//deep copy
	//align.globalGAP_diff = base_align.globalGAP_diff;//shallow copy
	//align.localGAP = base_align.localGAP;//shallow copy
	//auto localGAP = std::make_unique<std::vector<int32_t>>(sequence.size() + 1, 0);
	align.localGAP.assign(targetSequences.size(), std::vector<bool>(match.size() + 1, false));// + 1は番兵
	//align.localGAP.resize(targetSequences.size(), std::vector<int32_t>(LENGTH_MAX*2 + 1, false));
	//for (size_t mi = 0; mi < match.size(); ++mi)
	//{
	//	auto& npos = match[mi];

	//	if (s_aligni == npos.first) {
	//		//sGAP
	//		align.localGAP[align_sequence][mi] = true;
	//		for (size_t i = 0; i < targetSequences.size(); i++) {
	//			if (i != (size_t)align_sequence) {
	//				align.localGAP[i][mi] = base_align.localGAP[i][p_aligni];
	//			}
	//		}
	//		result_profile.push_back(profile[p_aligni]);
	//		result_profile.back().typenum[static_cast<int>(Type::GAP)] += 1;
	//	}
	//	else if (p_aligni == npos.second) {
	//		//pGAP
	//		for (size_t i = 0; i < targetSequences.size(); i++) {
	//			if (i != (size_t)align_sequence) {
	//				align.localGAP[i][mi] = true;
	//			}
	//		}
	//		result_profile.push_back(GapProfile);
	//		result_profile.back().typenum[static_cast<int>(sequence[s_aligni])] += 1;
	//	}
	//	else {
	//		for (size_t i = 0; i < targetSequences.size(); i++) {
	//			if (i != (size_t)align_sequence) {
	//				align.localGAP[i][mi] = base_align.localGAP[i][p_aligni];
	//			}
	//		}
	//		result_profile.push_back(profile[p_aligni]);
	//		result_profile.back().typenum[static_cast<int>(sequence[s_aligni])] += 1;
	//	}

	//	s_aligni = npos.first;
	//	p_aligni = npos.second;
	//}
	{
		s_aligni = 0;
		p_aligni = 0;
		for (size_t mi = 0; mi < match.size(); s_aligni = match[mi].first, p_aligni = match[mi].second, ++mi) {
			if (p_aligni == match[mi].second) {
				result_profile.push_back(GapProfile);
			}
			else {
				result_profile.push_back(profile[p_aligni]);
			}
			if (s_aligni == match[mi].first) {
				result_profile.back().typenum[static_cast<int>(Type::GAP)] += 1;
			}
			else {
				result_profile.back().typenum[static_cast<int>(sequence[s_aligni])] += 1;
			}
		}
	}
	{
		s_aligni = 0;
		for (size_t mi = 0; mi < match.size(); s_aligni = match[mi++].first) {
			if (s_aligni == match[mi].first) {
				align.localGAP[align_sequence][mi] = true;
			}
			else {
				//align.localGAP[align_sequence][mi] = false;
			}
		}
	}
	for (size_t i = 0; i < targetSequences.size(); i++) {
		if (i == (size_t)align_sequence) {
			continue;
		}
		p_aligni = 0;
		for (size_t mi = 0; mi < match.size(); p_aligni = match[mi++].second) {
			if (p_aligni == match[mi].second) {
				align.localGAP[i][mi] = true;
			}
			else {
				align.localGAP[i][mi] = base_align.localGAP[i][p_aligni];
			}
		}
	}
	//while (s_aligni > 0 && p_aligni > 0) {
	//	auto& dp_sprev = DP_MEMORY[s_aligni - 1];
	//	auto& dp = DP_MEMORY[s_aligni];
	//	++count;

	//	if (p_aligni > 0 && dp[p_aligni] == dp[p_aligni - 1] + calc_score(Type::GAP, profile[p_aligni - 1])) {
	//		//sGAP
	//		--p_aligni;
	//		if (profile[p_aligni].typenum == GapProfile.typenum) {
	//			--count;
	//			continue;
	//		}
	//		align.localGAP[align_sequence][profile_size - count] = true;
	//		result_profile.push_back(profile[p_aligni]);
	//		result_profile.back().typenum[static_cast<int>(Type::GAP)] += 1;
	//		continue;
	//	}
	//	assert(s_aligni > 0);
	//	if (dp[p_aligni] == dp_sprev[p_aligni] + calc_score(sequence[s_aligni - 1], GapProfile)) {
	//		//pGAP
	//		--s_aligni;
	//		for (size_t i = 0; i < targetSequences.size(); i++) {
	//			if (i != (size_t)align_sequence) {
	//				align.localGAP[i][profile_size - count] = true;
	//			}
	//		}
	//		result_profile.push_back(GapProfile);
	//		result_profile.back().typenum[static_cast<int>(sequence[s_aligni])] += 1;
	//		continue;
	//	}
	//	assert(p_aligni > 0);
	//	if (dp[p_aligni] == dp_sprev[p_aligni - 1] + calc_score(sequence[s_aligni - 1], profile[p_aligni - 1])) {
	//		//比較
	//		--p_aligni;
	//		--s_aligni;
	//		result_profile.push_back(profile[p_aligni]);
	//		result_profile.back().typenum[static_cast<int>(sequence[s_aligni])] += 1;
	//		continue;
	//	}
	//	assert(false);
	//}


	//return
	align.score = calc_profiles_score(result_profile);
	//align.globalGAP_diff[align_sequence] = align.globalGAP;
	//align.localGAP[align_sequence] = std::move(localGAP);

	//DEBUG
#ifdef _DEBUG
	{
		auto r_prof = align_to_array(align, align_sequence);
		if (r_prof.size() != result_profile.size()) {
			assert(r_prof.size() == result_profile.size());
		}
		r_prof = align_to_array(align, (align_sequence + 1) % targetSequences.size());
		if (r_prof.size() != result_profile.size()) {
			assert(r_prof.size() == result_profile.size());
		}
		r_prof = align_to_array(align, (align_sequence + 2) % targetSequences.size());
		if (r_prof.size() != result_profile.size()) {
			assert(r_prof.size() == result_profile.size());
		}
	}
#endif

	return { std::move(align), std::move(result_profile) };
}

void printAlign(const AlignData& align, const std::vector<Profile>& profile);
std::pair<AlignData, std::vector<Profile>> multiAlign(std::vector<std::vector<Type>> param_targetSequences)
{
	targetSequences = std::move(param_targetSequences);
	std::multimap<AlignData, std::vector<Profile>> queue;//探索候補
	std::multimap<AlignData, std::vector<Profile>> queue_next;
	std::mutex queue_next_mutex;
	{
		std::vector<Profile> start_profile;
		AlignData align;

		//初期化
		//align.globalGAP = std::make_unique<std::map<size_t, int32_t>>();
		//align.globalGAP_diff.resize(targetSequences.size(), align.globalGAP);
		//align.localGAP.resize(targetSequences.size());
		auto max_length_seq = std::max_element(all_range(targetSequences), [](auto& l, auto& r) {return l.size() < r.size(); });
		start_profile.resize(max_length_seq->size());
		align.localGAP.assign(targetSequences.size(), std::vector<bool>(start_profile.size() + 1, false));
		for (size_t si = 0; si < targetSequences.size(); ++si)
		{

			////後ろにGAP
			//auto& seq = targetSequences[si];
			//size_t i = 0;
			//for (; i < seq.size(); ++i) {
			//	start_profile[i].typenum[static_cast<int>(seq[i])] += 1;
			//}
			////align.localGAP[si] = std::make_unique<std::vector<int32_t>>(seq.size() + 1);
			////(*align.localGAP[si])[seq.size()] = start_profile.size() - i;
			//for (; i < start_profile.size(); ++i) {
			//	start_profile[i].typenum[static_cast<int>(Type::GAP)] += 1;
			//	align.localGAP[si][i] = true;
			//}

			//最大長さにアライン
			if (static_cast<int>(si) != max_length_seq - targetSequences.begin()) {
				SequenceFixedSequenceAlign(si, static_cast<int>(max_length_seq - targetSequences.begin()), align);
			}
			auto&& aligned_seq = align_to_array(align, si);
			for (size_t i = 0; i < max_length_seq->size(); i++)
			{
				start_profile[i].typenum[static_cast<int>(aligned_seq[i])] += 1;
			}
		}
		align.score = calc_profiles_score(start_profile);

		queue.emplace(std::move(align), std::move(start_profile));
	}

	//ビームサーチ
#ifdef _DEBUG
	static constexpr int32_t width = 3;
	static constexpr int32_t random_num = 3;
	static constexpr int32_t LOOP_NUM = 3;
#else
	static constexpr int32_t width = 25;
	static constexpr int32_t random_num = 4;
	static constexpr int32_t LOOP_NUM = 1000;
#endif
	static int_fast32_t DP_MEMORY[width][random_num][LENGTH_MAX + 1][2 * LENGTH_MAX + 1];
	static int_fast32_t last10_score = -10 * LENGTH_SUM_MAX;

	for (size_t COUNTCOUNT = 1; COUNTCOUNT <= LOOP_NUM; COUNTCOUNT++)
	{
		if (COUNTCOUNT % 10 == 0) {
			std::cout << COUNTCOUNT << std::endl;
			auto newscore = queue.rbegin()->first.score;
			if (last10_score == newscore) {
				std::cout << "NOT UPDATED" << '\n';
				break;
			}
			last10_score = newscore;
			if (COUNTCOUNT % 100 == 0) {
				std::cout << "========================" << '\n';
				printAlign(queue.rbegin()->first, queue.rbegin()->second);
				std::cout << "========================" << std::endl;
			}
		}

		std::for_each(std::execution::par, make_indexed_iterator(queue.begin()), make_indexed_iterator(queue.end()),
			[&queue_next, &queue_next_mutex](auto base_with_index) {
				auto& base = base_with_index.second;
				auto& DP_MEMORY_LOCAL = DP_MEMORY[base_with_index.first];
				//アライメントをするシーケンスの番号をランダムに重複なく選択
				std::vector<int> selected_seq(targetSequences.size());
				std::iota(all_range(selected_seq), 0);
				std::shuffle(all_range(selected_seq), std::mt19937{ std::random_device{}() });
				selected_seq.resize(std::min<size_t>(random_num, targetSequences.size()));
				selected_seq.shrink_to_fit();

				std::for_each(std::execution::par, make_indexed_iterator(selected_seq.begin()), make_indexed_iterator(selected_seq.end()),
					[&base, &DP_MEMORY_LOCAL, &queue_next, &queue_next_mutex](auto r) {

						std::vector<Profile> profile = base.second;//copy
						{
							//profileから選択したものを除外
							const auto& r_prof = align_to_array(base.first, r.second);
							assert(r_prof.size() == profile.size());
							for (size_t i = 0; i < r_prof.size(); i++) {
								profile[i].typenum[static_cast<int>(r_prof[i])] -= 1;
							}
						}
						auto&& nextAlign = SequenceProfileAlign(r.second, profile, base.first, DP_MEMORY_LOCAL[r.first]);
						assert(base.first.score - std::abs(base.first.score) * 0.001 < nextAlign.first.score);
						assert(nextAlign.first.score == calc_profiles_score(nextAlign.second));
						//追加
						std::lock_guard<std::mutex> lock(queue_next_mutex);
						queue_next.emplace(nextAlign);
					});
			});

		//幅制限
		{
			auto iter = queue_next.begin();
			for (size_t i = 0; i < width && iter != queue_next.end(); i++) {
				++iter;
			}
			queue_next.erase(iter, queue_next.end());
		}
		//queue = std::move(queue_next);
		//queue_next.clear();
		queue = queue_next;//copy 最適値が消えるのを防ぐため
	}

	return std::move(*queue.rbegin());
}


constexpr Type char_to_Type(char t) {
	if ('a' <= t && t <= 'z') {
		t = t - 'a' + 'A';
	}
	switch (t)
	{
	case 'A': return Type::A;
	case 'R': return Type::R;
	case 'N': return Type::N;
	case 'D': return Type::D;
	case 'C': return Type::C;
	case 'Q': return Type::Q;
	case 'E': return Type::E;
	case 'G': return Type::G;
	case 'H': return Type::H;
	case 'I': return Type::I;
	case 'L': return Type::L;
	case 'K': return Type::K;
	case 'M': return Type::M;
	case 'F': return Type::F;
	case 'P': return Type::P;
	case 'S': return Type::S;
	case 'T': return Type::T;
	case 'W': return Type::W;
	case 'Y': return Type::Y;
	case 'V': return Type::V;
	case 'B': return Type::B;
	case 'Z': return Type::Z;
	case 'X': return Type::X;
	case '-': return Type::GAP;
	default:
		return Type::GAP;
	}
}
constexpr char Type_to_char(Type t) {
	constexpr char table[] = { 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','-' };
	return table[static_cast<int>(t)];
}
std::vector<Type> str_to_Type(std::string_view str) {
	std::vector<Type> res;
	res.reserve(str.size());
	for (auto& c : str) {
		res.push_back(char_to_Type(c));
	}
	return res;
}

void printAlign(const AlignData& align, const std::vector<Profile>& profile)
{
	std::cout << "score = " << align.score << '\n';
	std::string res;
	{
		res.clear();
		std::string score_str;
		for (auto& pr : profile) {
			//平均スコア3.5以上
			auto score = calc_profile_score(pr);
			if (score >=
					static_cast<int_fast32_t>(7 * align.localGAP.size() * (align.localGAP.size() - 1) / 2 / 2)) {
				res += '*';
			}
			else {
				res += ' ';
			}

			if (score < 0) {
				score_str += '-';
			}
			else if (score >= 10) {
				score_str += '+';
			}
			else {
				score_str += (char)('0' + score);
			}
		}
		std::cout << res << '\n';
	}
	for (size_t si = 0; si < targetSequences.size(); si++)
	{
		res.clear();
		for (auto& t : align_to_array(align, si)) {
			res += Type_to_char(t);
		}
		std::cout << res << '\n';
	}
	std::cout << std::flush;
}

void testD()
{
	const char* const testcases_str[] = {
		"mkpiiavnykayypysfgenalriardakrvweetgvevilappfteiyrvlkevegsgvkvfaqhavmkdmalnalkavs",
		"maprkffvggnwkmngdkkslgelihtlngaklsadtevvcgapsiyldfarqkldakigvaaqncyckelasqhdvdgflvggaslkpefvdiinakh",
		"mrryliagnwkmntsletgtalasgladhvrgrdlpvdvlvcppfpylaavkatageagisvgaqncyfeasgaftgevsvdmlkdigcds",
	};

	std::vector<std::vector<Type>> target;
	target.reserve(sizeof(testcases_str) / sizeof(testcases_str[0]));
	for (auto& s : testcases_str) {
		target.emplace_back(str_to_Type(s));
	}

	auto&& res = multiAlign(std::move(target));
	printAlign(res.first, res.second);
}
void test1()
{
	//このケースは友人に頂きました
	const char* const testcases_str[] = {
		//"mkpiiavnykayypysfgenalriardakrvweetgvevilappfteiyrvlkevegsgvkvfaqhavmkdmalnalkavs",
		//"maprkffvggnwkmngdkkslgelihtlngaklsadtevvcgapsiyldfarqkldakigvaaqncyckelasqhdvdgflvggaslkpefvdiinakh",
		//"mrryliagnwkmntsletgtalasgladhvrgrdlpvdvlvcppfpylaavkatageagisvgaqncyfeasgaftgevsvdmlkdigcds",

		"mkpiiavnykayypysfgenalriardakrvweetgvevilappfteiyrvlkevegsgvkvfaqhadpvepgavtgyipveglkeagvhgvilnhsehrlkiadinaliikarrlglktlacadvpetgaaiallkpdmiaveppeligtgvsvskakpevitnsvsmirsvnkealiltgagittgedvyqavklgtigvlvasgivkakdpysvmkdmalnalkavs",
		"maprkffvggnwkmngdkkslgelihtlngaklsadtevvcgapsiyldfarqkldakigvaaqncykvpkgaftgeispamikdigaawvilghserrhvfgesdeligqkvahalaeglgviacigekldereagitekvvfeqtkaiadnvkdwskvvlayepvwaigtgktatpqqaqevheklrgwlkshvsdavaqstriiyggsvtggnckelasqhdvdgflvggaslkpefvdiinakh",
		"maatsltappsfsglrrispkldaaavsshqsffhrvnsstrlvsssssshrsprgvvamagsgkngtkdsiaklisdlnsatleadvdvvvsppfvyidqvkssltdridisgqnswvgkggaftgeisveqlkdlgckwvilghserrhvigekdefigkkaayalseglgviacigekleereagktfdvcfaqlkafadavpswdnivvayepvwaigtgkvaspqqaqevhvavrgwlkknvseevasktriiyggsvnggnsaelakeedidgflvggaslkgpefativnsvtskkvaa",
		"markffvggnwkcngtaeevkkivntlneaqvpsqdvvevvvsppyvflplvkstlrsdffvaaqncwvkkggaftgevsaemlvnldipwvilghserrailnessefvgdkvayalaqglkviacvgetleereagstmdvvaaqtkaiadrvtnwsnvviayepvwaigtgkvaspaqaqevhdelrkwlaknvsadvaattriiyggsvnggnckelggqadvdgflvggaslkpefidiikaaevkksa",
		"mrqiiiagnwkmhktqtesleflqgflshledtpeeretvlcvpftclnfmsknlhgsrvklgaqnihwadqgaftgeisgemlkefginyvivghserrqyfgetdetvnarllaaqkhgltpilcvgeskaqrdageteavisaqiekdlvnvdqnnlviayepiwaigtgdtceaaeanrviglirsqltnknvtiqyggsvnpknvdeimaqpeidgalvggasldpesfarlvnyq",
		"mhktqaesleflqsflpqlentaedrevilcapytalgvmsknlhgtrvrigsqnvhweesgaftgeiapsmlteigvtyavvghserrqyfgetdetvnfraraaqkaeltpilcvgeskeqrdagqtetvikeqlkadlvgvdlsqlviayepiwaigtgdtceaeeanrvigmirselsssdvpiqyggsvkpanideimaqpeidgalvggasldpvgfarivnyeat",
		"mkrqiviagnwkmhktnseamqlanqvriktmditktqivicppftalapvyevigdsrihlgaqnmfwekegaftgeisagmikstgadyviighserrqyfgesdetvnkkvkaalenglkpivcvgetleereanitlkvvsrqirgafadlsaeqmkkvivayepvwaigtgktatpeqaqqvhqeirqlltemfgseigekmviqyggsvkpanaesllsqpdidgalvggaclkadsfseiihiaeklq",
		"mktrqqivagnwkmnknygegrelameiverlkpsntqvvlcapyihlqlvkniikdvaslylgaqnchqedkgaytgeisvdmlksvgvsyvilghserreyfgesdellakktdkvlaagllpifccgesldirdagthvahvqaqikaglfhlspeefqkvviayepiwaigtgrtaspeqaqdmhaairalltdqygaeiadattilyggsvnggnaavlfsqpdvdgglvggaslkaeefitiveatkk",
		"mvywvgtswkmnktlaeamdfaailagfvpgfddriqpfvippftavrqvkqalsstrvkvgaqnmhwadagawtgeispvmltdcgldlvelghserrehfgetdrtvglktaaavkhgliplicvgetlaeresgeadavlakqvegalqffeeevkgatilfayepvwaigdkgipassdyadkqqglikavagsllpsvpsvlyggsvnpgnaaeligqpnvdglfigrsawqaqgyidilgrasaai",
		"mrryliagnwkmntsletgtalasgladhvrgrdlpvdvlvcppfpylaavkatageagisvgaqncyfeasgaftgevsvdmlkdigcdsvilghserrhvikecddminkktkaaiegglqvvlcvgelleereadkteavldeqmagglkdisaeqmtnvviayepvwaigtgktaspeqaeqahahlrkwladrytsevaeqtrilyggsvkpanakellgqqnvdgalvggasltvdnfgpiidagvelsa",
		"mpeekpviminfktynesygfrahdiaeaaetvaeesgieivicpgfmdihpmsnhyrlpvfaqhidgispgahtghilaeavraagatgtlinhserrltladisaavdaakranlktvvctnntatsgaaaalspdyvaieppeligsgisvatadpeiiensvnavksvnkdvkvlagagissgscvkravelgsdgvllasgvvkaedpavvlrdlvski",
		"mgsplivvnfktylegtgersvdiaracrdvaedsgvdiavapqmcdiyrvasmvdipvysqhvdgigagsftghafapaikeagasgtlinhsenrltladieaaiqaskavglktivctnniptsaaaaalspdyvaveppeligsgipvseadpdvvkgsveavmnidsgvsvlcgagiskgkdlkaaldlgskgvllasgivksedprsamedlisli",
		"mrkkivagnwkmnldytegltlfsevinmikdevtgsqqavvcspfihlhslvqlgkdynkvsvgaqnahqaeagaytgeissrmiksvgaeyviighserrqyfgetndllakktdavlknqltpifcigetlqeretekhfeviksqllegvfhldetafaklviayepvwaigtgvtasaeqaqeihafiraeiaqkysqqvadditilyggscnpknaaelfakgdidggliggaslksrdfvdilkvfn",
		"mdkleakesaclslsssriggmrkkliagnwkmnqtpsqavvladalkktvsgkeaaeivvcppytalipvrdalkgssihlgaqdlhwedqgaftgkisadmlldagcthviighseqrtyfhetdatvnkklikalagglvpifcigetleerdggrafdvvkkqleggfagmkdaghtvlayepvwaigtgrnatpeqaqemhafirktiaslfsaavadgmrilyggsmkpdnaagllaqpdidggliggaalkadsfygivkaag",
		"mtpaapgapavqrrplfagnwkmhtlpaeaarlaaavregldgwgggdpaaaggtvtpgvagrgagtqpaegpagppaagaaevvlcppftslaaaaealagsaialgaqdlawgdfgaftgevsapmlrelgcryvivghserrqllgetdalilrkleaalagglvpilcvgedaaqrregrtaavvlgqaahalagldgeqaarvviayepvwaigsgtpatpadaqavaaalrglierlhgpavaaavrilyggsvkpdniggfmaqpdidgalvggasldgagfarlvrqgvaaraaapggaaageers",
		"matsktvgrvplmagnwkmnldhlqathliqkldwtlrdakhdydgvevavlppftdlrsvqtlvegdrlhlrygaqdlsphasgaytgdisgaflkklgctyvvvghserreghhetddvvaakvqaayrhgltpilccgeglevrkegsqvehvvaqlraaldgvtreqaasiviayepiwaigtgevatpddaqevcaairtllaelysgdladgvrilyggsvkaanvaaimaqedvdgalvggasidpaefasicryrdhltag",
		"mrtkmiagnwkmhhrpqearafveelgrvlwarnelygplkegvaeavlfptalslaavqdalgdlpvslgaqnahwedhgaftgeigapmladfgcayilighserrhlfhetevelarklravlstsarclfcvgelleereagkthqvlerqllgalegvtipdltdrfavayepvwaigtgktasdgdaeegcgyirhlvadrygqetaqhlqvlyggsvkpgntaglmvqgdidgllvggaslevpsfvgileaalgilrp",
	};

	std::vector<std::vector<Type>> target;
	target.reserve(sizeof(testcases_str) / sizeof(testcases_str[0]));
	for (auto& s : testcases_str) {
		target.emplace_back(str_to_Type(s));
	}

	auto&& res = multiAlign(std::move(target));
	printAlign(res.first, res.second);
}
void test2()
{
	//このケースは友人に頂きました
	const char* const testcases_str[] = {
"mkpiiavnykayypysfgenalriardakrvweetgvevilappfteiyrvlkevegsgvkvfaqhadpvepgavtgyipveglkeagvhgvilnhsehrlkiadinaliikarrlglktlacadvpetgaaiallkpdmiaveppeligtgvsvskakpevitnsvsmirsvnkealiltgagittgedvyqavklgtigvlvasgivkakdpysvmkdmalnalkavs",
"maprkffvggnwkmngdkkslgelihtlngaklsadtevvcgapsiyldfarqkldakigvaaqncykvpkgaftgeispamikdigaawvilghserrhvfgesdeligqkvahalaeglgviacigekldereagitekvvfeqtkaiadnvkdwskvvlayepvwaigtgktatpqqaqevheklrgwlkshvsdavaqstriiyggsvtggnckelasqhdvdgflvggaslkpefvdiinakh",
"markffvggnwkcngtaeevkkivntlneaqvpsqdvvevvvsppyvflplvkstlrsdffvaaqncwvkkggaftgevsaemlvnldipwvilghserrailnessefvgdkvayalaqglkviacvgetleereagstmdvvaaqtkaiadrvtnwsnvviayepvwaigtgkvaspaqaqevhdelrkwlaknvsadvaattriiyggsvnggnckelggqadvdgflvggaslkpefidiikaaevkksa",
"mrqiiiagnwkmhktqtesleflqgflshledtpeeretvlcvpftclnfmsknlhgsrvklgaqnihwadqgaftgeisgemlkefginyvivghserrqyfgetdetvnarllaaqkhgltpilcvgeskaqrdageteavisaqiekdlvnvdqnnlviayepiwaigtgdtceaaeanrviglirsqltnknvtiqyggsvnpknvdeimaqpeidgalvggasldpesfarlvnyq",
"mhktqaesleflqsflpqlentaedrevilcapytalgvmsknlhgtrvrigsqnvhweesgaftgeiapsmlteigvtyavvghserrqyfgetdetvnfraraaqkaeltpilcvgeskeqrdagqtetvikeqlkadlvgvdlsqlviayepiwaigtgdtceaeeanrvigmirselsssdvpiqyggsvkpanideimaqpeidgalvggasldpvgfarivnyeat",
"mkrqiviagnwkmhktnseamqlanqvriktmditktqivicppftalapvyevigdsrihlgaqnmfwekegaftgeisagmikstgadyviighserrqyfgesdetvnkkvkaalenglkpivcvgetleereanitlkvvsrqirgafadlsaeqmkkvivayepvwaigtgktatpeqaqqvhqeirqlltemfgseigekmviqyggsvkpanaesllsqpdidgalvggaclkadsfseiihiaeklq",
"mktrqqivagnwkmnknygegrelameiverlkpsntqvvlcapyihlqlvkniikdvaslylgaqnchqedkgaytgeisvdmlksvgvsyvilghserreyfgesdellakktdkvlaagllpifccgesldirdagthvahvqaqikaglfhlspeefqkvviayepiwaigtgrtaspeqaqdmhaairalltdqygaeiadattilyggsvnggnaavlfsqpdvdgglvggaslkaeefitiveatkk",
"mvywvgtswkmnktlaeamdfaailagfvpgfddriqpfvippftavrqvkqalsstrvkvgaqnmhwadagawtgeispvmltdcgldlvelghserrehfgetdrtvglktaaavkhgliplicvgetlaeresgeadavlakqvegalqffeeevkgatilfayepvwaigdkgipassdyadkqqglikavagsllpsvpsvlyggsvnpgnaaeligqpnvdglfigrsawqaqgyidilgrasaai",
"mrryliagnwkmntsletgtalasgladhvrgrdlpvdvlvcppfpylaavkatageagisvgaqncyfeasgaftgevsvdmlkdigcdsvilghserrhvikecddminkktkaaiegglqvvlcvgelleereadkteavldeqmagglkdisaeqmtnvviayepvwaigtgktaspeqaeqahahlrkwladrytsevaeqtrilyggsvkpanakellgqqnvdgalvggasltvdnfgpiidagvelsa",
"mpeekpviminfktynesygfrahdiaeaaetvaeesgieivicpgfmdihpmsnhyrlpvfaqhidgispgahtghilaeavraagatgtlinhserrltladisaavdaakranlktvvctnntatsgaaaalspdyvaieppeligsgisvatadpeiiensvnavksvnkdvkvlagagissgscvkravelgsdgvllasgvvkaedpavvlrdlvski",
"mgsplivvnfktylegtgersvdiaracrdvaedsgvdiavapqmcdiyrvasmvdipvysqhvdgigagsftghafapaikeagasgtlinhsenrltladieaaiqaskavglktivctnniptsaaaaalspdyvaveppeligsgipvseadpdvvkgsveavmnidsgvsvlcgagiskgkdlkaaldlgskgvllasgivksedprsamedlisli",
"mrkkivagnwkmnldytegltlfsevinmikdevtgsqqavvcspfihlhslvqlgkdynkvsvgaqnahqaeagaytgeissrmiksvgaeyviighserrqyfgetndllakktdavlknqltpifcigetlqeretekhfeviksqllegvfhldetafaklviayepvwaigtgvtasaeqaqeihafiraeiaqkysqqvadditilyggscnpknaaelfakgdidggliggaslksrdfvdilkvfn",
"matsktvgrvplmagnwkmnldhlqathliqkldwtlrdakhdydgvevavlppftdlrsvqtlvegdrlhlrygaqdlsphasgaytgdisgaflkklgctyvvvghserreghhetddvvaakvqaayrhgltpilccgeglevrkegsqvehvvaqlraaldgvtreqaasiviayepiwaigtgevatpddaqevcaairtllaelysgdladgvrilyggsvkaanvaaimaqedvdgalvggasidpaefasicryrdhltag",
"mrtkmiagnwkmhhrpqearafveelgrvlwarnelygplkegvaeavlfptalslaavqdalgdlpvslgaqnahwedhgaftgeigapmladfgcayilighserrhlfhetevelarklravlstsarclfcvgelleereagkthqvlerqllgalegvtipdltdrfavayepvwaigtgktasdgdaeegcgyirhlvadrygqetaqhlqvlyggsvkpgntaglmvqgdidgllvggaslevpsfvgileaalgilrp",
	};

	std::vector<std::vector<Type>> target;
	target.reserve(sizeof(testcases_str) / sizeof(testcases_str[0]));
	for (auto& s : testcases_str) {
		target.emplace_back(str_to_Type(s));
	}

	auto&& res = multiAlign(std::move(target));
	printAlign(res.first, res.second);
}
int main()
{
#ifdef _DEBUG
	testD();
	return 0;
#endif
	test1();
	//test2();

	return 0;
}
#endif
