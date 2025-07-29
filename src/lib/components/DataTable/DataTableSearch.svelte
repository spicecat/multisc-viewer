<script lang="ts">
	import IconButton from '@smui/icon-button';
	import TextField from '@smui/textfield';
	import Fuse, { type IFuseOptions } from 'fuse.js';
	import DataTable from './DataTable.svelte';

	let {
		label = '',
		data,
		columns,
		isLoading,
		selected = $bindable(),
		searchOptions = {}
	}: {
		label: string;
		data: Data[];
		columns: Column[];
		isLoading?: boolean;
		selected?: string | string[];
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

<div style="display: flex;flex-direction:column;">
	<div>
		<TextField bind:value={query} {label}>
			{#snippet leadingIcon()}
				<IconButton class="material-icons">search</IconButton>
			{/snippet}
		</TextField>
	</div>
	<DataTable data={filteredItems} {columns} {isLoading} bind:selected />
</div>
