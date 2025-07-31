<script lang="ts">
	import { makeTitle } from '$lib/utils/utils';
	import Checkbox from '@smui/checkbox';
	import Chip, { Set, Text } from '@smui/chips';
	import DataTable, { Body, Cell, Head, Label, Pagination, Row } from '@smui/data-table';
	import IconButton from '@smui/icon-button';
	import LinearProgress from '@smui/linear-progress';
	import Radio from '@smui/radio';
	import Select, { Option } from '@smui/select';
	import type { Data } from '$lib/types/data';
	import type { Columns } from '$lib/types/data-table';

	let {
		data,
		columns,
		isLoading = false,
		selected = $bindable()
	}: {
		data: Data[];
		columns: Columns;
		isLoading?: boolean;
		selected?: string | string[];
	} = $props();

	let sort = $state('');
	let sortDirection: 'ascending' | 'descending' = $state('ascending');

	const id = $derived(columns[0]?.key);

	const items = $derived(
		[...data]
			.map((item) =>
				item && typeof item === 'object' && !Array.isArray(item) ? item : { '': item }
			)
			.sort((a, b) => {
				const [aVal, bVal] = [a[sort], b[sort]][
					sortDirection === 'ascending' ? 'slice' : 'reverse'
				]();
				if (typeof aVal === 'string' && typeof bVal === 'string') return aVal.localeCompare(bVal);
				return aVal > bVal ? 1 : aVal < bVal ? -1 : 0;
			})
	);

	let perPage = $state(10);
	let currentPage = $state(0);
	const start = $derived(currentPage * perPage);
	const end = $derived(Math.min(start + perPage, items.length));
	const slice = $derived(items.slice(start, end));
	const lastPage = $derived(Math.max(Math.ceil(items.length / perPage) - 1, 0));
	$effect(() => {
		if (currentPage > lastPage) currentPage = lastPage;
	});
</script>

<DataTable sortable stickyHeader bind:sort bind:sortDirection table$aria-label="Data table">
	<Head>
		<Row>
			{#if selected !== undefined}
				<Cell checkbox>
					{#if Array.isArray(selected)}
						<Checkbox />
					{/if}
				</Cell>
			{/if}
			{#each columns as { key, label, sortable = true }}
				<Cell columnId={key} numeric={typeof items[0]?.[key] === 'number'} {sortable}>
					<Label>{label}</Label>
					{#if sortable}
						<IconButton class="material-icons">arrow_upward</IconButton>
					{/if}
				</Cell>
			{/each}
		</Row>
	</Head>
	{#snippet progress()}
		<LinearProgress indeterminate closed={!isLoading} aria-label="Data is being loaded..." />
	{/snippet}
	<Body style="white-space: normal;">
		{#each slice as item (item[id])}
			<Row>
				{#if selected !== undefined}
					<Cell checkbox>
						{#if Array.isArray(selected)}
							<Checkbox bind:group={selected} value={item[id]} />
						{:else}
							<Radio bind:group={selected} value={item[id]} />
						{/if}
					</Cell>
				{/if}
				{#each columns as { key, url }}
					<Cell numeric={typeof item[key] === 'number'}>
						{#if url}
							<a href={item[key].href}>{item[key].name}</a>
						{:else if Array.isArray(item[key])}
							<Set chips={item[key]}>
								{#snippet chip(chip)}
									<Chip {chip}>
										<Text>{chip}</Text>
									</Chip>
								{/snippet}
							</Set>
						{:else}
							{key === 'name' ? makeTitle(item[key]) : item[key]}
						{/if}
					</Cell>
				{/each}
			</Row>
		{/each}
	</Body>
	{#snippet paginate()}
		<Pagination>
			{#snippet rowsPerPage()}
				<Label>Rows Per Page</Label>
				<Select variant="outlined" bind:value={perPage} noLabel>
					<Option value={10}>10</Option>
					<Option value={25}>25</Option>
					<Option value={100}>100</Option>
				</Select>
			{/snippet}
			{#snippet total()}
				{start + 1}-{end} of {items.length}
			{/snippet}

			<IconButton
				class="material-icons"
				action="first-page"
				title="First page"
				onclick={() => (currentPage = 0)}
				{...{ disabled: currentPage === 0 }}
			>
				first_page
			</IconButton>
			<IconButton
				class="material-icons"
				action="prev-page"
				title="Prev page"
				onclick={() => currentPage--}
				{...{ disabled: currentPage === 0 }}
			>
				chevron_left
			</IconButton>
			<IconButton
				class="material-icons"
				action="next-page"
				title="Next page"
				onclick={() => currentPage++}
				{...{ disabled: currentPage === lastPage }}
			>
				chevron_right
			</IconButton>
			<IconButton
				class="material-icons"
				action="last-page"
				title="Last page"
				onclick={() => (currentPage = lastPage)}
				{...{ disabled: currentPage === lastPage }}
			>
				last_page
			</IconButton>
		</Pagination>
	{/snippet}
</DataTable>
