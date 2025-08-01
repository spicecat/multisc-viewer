<script lang="ts">
	import type { Columns, Select, Data } from '$lib/types/data-table';
	import { Search } from '@lucide/svelte';
	import Fuse, { type IFuseOptions } from 'fuse.js';
	import DataTable from './DataTable.svelte';

	let {
		label = '',
		data,
		columns,
		select,
		searchOptions = {}
	}: {
		label: string;
		data: Data[];
		columns: Columns;
		select: Select;
		searchOptions?: IFuseOptions<Data>;
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

	const filteredItems = $derived(query ? fuse.search(query).map(({ item }) => item) : data);
</script>

<div class="input-group grid-cols-[auto_1fr_auto]">
	<div class="ig-cell preset-tonal">
		<Search size={16} />
	</div>
	<input class="ig-input" type="text" bind:value={query} placeholder={`Search ${label}...`} />
</div>
<DataTable data={filteredItems} {columns} {select} />
