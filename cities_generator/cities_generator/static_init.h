#pragma once

EventHandler* EventHandler::handler = nullptr;
ShaderLibrary* ShaderLibrary::library = nullptr;
unsigned int Entity::currentId = 0;
std::vector<Entity*> Entity::entities = std::vector<Entity*>{};
std::unordered_map<unsigned int, Entity*> Entity::entities_map = std::unordered_map<unsigned int, Entity*>{};
std::array<std::unordered_set<unsigned>, COMPONENT_COUNTER> Component::componentOwnersArray{};
std::map<std::string, SleepInterval*> SleepInterval::intervals{};
std::unordered_map<std::string,std::unordered_set<unsigned>> TagComponent::tagOwnersMap;
// std::vector<Archetype> Archetype::archetypes = std::vector<Archetype>{};
Archetype* Archetype::basicNode = nullptr;
std::unordered_map<unsigned int, Record> Entity::entitiesArchetype = std::unordered_map<unsigned int, Record>{};
bool Archetype::hasUnhandledArchetypes = false;
std::vector<System*> System::allSystems =  std::vector<System*>{};