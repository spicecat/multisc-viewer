<script lang="ts">
	import Radio from '@smui/radio';
	import FormField from '@smui/form-field';

	let { groupBy = $bindable() }: { groupBy: Grouping } = $props();

	const groupOptions = [
		{ value: 'Genotype', label: 'Group by Genotype' },
		{ value: 'CellType', label: 'Group by Cell Type' }
	];
</script>

<div>
	{#each groupOptions as { value, label }}
		<FormField>
			<Radio bind:group={groupBy} {value} />
			<span>{label}</span>
		</FormField>
	{/each}
</div>
