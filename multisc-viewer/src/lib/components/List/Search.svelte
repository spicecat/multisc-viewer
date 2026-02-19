<script lang="ts" generics="T">
import { TagsInput } from "@skeletonlabs/skeleton-svelte";
import Fuse, { type IFuseOptions } from "fuse.js";

let {
	name,
	data,
	items = $bindable(),
	tags = $bindable(),
	searchOptions = {},
}: {
	name: string;
	data: T[];
	items: T[];
	tags: string[];
	searchOptions?: IFuseOptions<T>;
} = $props();

const defaultSearchOptions = $derived({
	threshold: 0.0,
	ignoreLocation: true,
	useExtendedSearch: true,
});

const fuse = $derived(
	new Fuse(data, {
		...defaultSearchOptions,
		...searchOptions,
	}),
);

$effect(() => {
	items = tags.length
		? fuse
				.search({
					$or: tags.map((tag) => {
						const [val, path] = tag.toLowerCase().split("=", 2).toReversed();
						return {
							$or:
								searchOptions.keys
									?.filter((key) => {
										const keyStr = key.toString().toLowerCase();
										return (
											!path || keyStr === path || keyStr.startsWith(`${path}.`)
										);
									})
									.map((path) => ({
										$path: path.toString(),
										$val: val,
									})) ?? [],
						};
					}),
				})
				.map(({ item }) => item)
		: data;
});
</script>

<TagsInput value={tags} addOnPaste name="tags" onValueChange={(e) => (tags = e.value)}>
	<TagsInput.Label>Search {name}</TagsInput.Label>
	<TagsInput.Control>
		<TagsInput.Context>
			{#snippet children(tagsInput)}
				{#each tagsInput().value as value, index (index)}
					<TagsInput.Item {value} {index}>
						<TagsInput.ItemPreview class="preset-tonal-primary">
							<TagsInput.ItemText>{value}</TagsInput.ItemText>
							<TagsInput.ItemDeleteTrigger />
						</TagsInput.ItemPreview>
						<TagsInput.ItemInput />
					</TagsInput.Item>
				{/each}
			{/snippet}
		</TagsInput.Context>
		<TagsInput.Input placeholder="Add a tag..." />
	</TagsInput.Control>
	<TagsInput.ClearTrigger>Clear All</TagsInput.ClearTrigger>
	<TagsInput.HiddenInput />
</TagsInput>
