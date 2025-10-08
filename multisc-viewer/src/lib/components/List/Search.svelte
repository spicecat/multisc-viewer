<script lang="ts" generics="T">
	import { Search } from '@lucide/svelte';
	import { TagsInput } from '@skeletonlabs/skeleton-svelte';
	import Fuse, { type IFuseOptions } from 'fuse.js';

	let {
		name,
		data,
		items = $bindable(),
		tags = $bindable(),
		searchOptions = {}
	}: {
		name: string;
		data: T[];
		items: T[];
		tags: string[];
		searchOptions?: IFuseOptions<T>;
	} = $props();

	let defaultSearchOptions = $derived({
		threshold: 0.0,
		ignoreLocation: true,
		useExtendedSearch: true
	});

	let fuse = $derived(
		new Fuse(data, {
			...defaultSearchOptions,
			...searchOptions
		})
	);

	$effect(() => {
		items = tags.length
			? fuse
					.search({
						$and: tags.map((tag) => {
							const [val, path] = tag.toLowerCase().split('=', 2).toReversed();
							return {
								$or:
									searchOptions.keys
										?.filter((key) => {
											const keyStr = key.toString().toLowerCase();
											return !path || keyStr === path || keyStr.startsWith(path + '.');
										})
										.map((path) => ({
											$path: path.toString(),
											$val: val
										})) ?? []
							};
						})
					})
					.map(({ item }) => item)
			: data;
	});
</script>

<span class="input-group flex w-full">
	<div class="ig-cell preset-tonal">
		<Search size={16} />
	</div>
	<TagsInput
		name="tags"
		value={tags}
		onValueChange={(e) => (tags = e.value)}
		placeholder={`Search ${name}...`}
	/>
</span>
