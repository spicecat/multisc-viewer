<script lang="ts">
	import type { Columns, Data, Select } from '$lib/types/data-table';
	import { Search } from '@lucide/svelte';
	import Fuse, { type IFuseOptions } from 'fuse.js';
	import type { Snippet } from 'svelte';
	import DataTable from './DataTable.svelte';

	let {
		name = '',
		data,
		columns,
		select,
		searchOptions = {},
		children
	}: {
		name: string;
		data: Data[];
		columns: Columns;
		select?: Select;
		searchOptions?: IFuseOptions<Data>;
		children?: Snippet;
	} = $props();

	const defaultSearchOptions = $derived({
		keys: columns.map(({ key }) => key),
		threshold: 0.0,
		ignoreLocation: true,
		useExtendedSearch: true
	});

	const fuse = $derived(
		new Fuse(data, {
			...defaultSearchOptions,
			...searchOptions
		})
	);

	let query = $state('');

	const items = $derived(query ? fuse.search(query).map(({ item }) => item) : data);
</script>

<div class="space-y-2">
	<div class="input-group grid-cols-[auto_1fr_auto]">
		<div class="ig-cell preset-tonal">
			<Search size={16} />
		</div>
		<input class="ig-input" type="text" bind:value={query} placeholder={`Search ${name}...`} />
	</div>
	{@render children?.()}
	<DataTable {name} data={items} {columns} {select} />
</div>
