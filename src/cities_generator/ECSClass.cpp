#include "cities_generator/ECSClass.h"
#include <algorithm>
#include <iostream>

const bool IS_ARCHETYPE_DEBUG = false;

Entity::Entity()
{
    id = Entity::currentId++;
    Entity::entities_map.emplace(id, this);

    Entity::entities.push_back(this);

    Archetype::AddEmptyEntity(id);
}

Entity::~Entity()
{
    // for (size_t i = 0; i < components.size(); i++)
    // {
    //     delete components[i];
    // }
    Entity::entities.erase(std::find(Entity::entities.begin(),Entity::entities.end(), this));
    Entity::entities_map.erase(id);

    Record* cur = &Entity::entitiesArchetype[id];
    for (int i = 0; i < cur->archetype->components.size(); i++)
    {
        delete cur->archetype->components[i][cur->index];
        cur->archetype->components[i].erase(cur->archetype->components[i].begin() + cur->index);
    }
    cur->archetype->entitiesId.erase(cur->archetype->entitiesId.begin() + cur->index);
    for (int i = cur->index; i < cur->archetype->entitiesId.size(); i++)
    {
        Entity::entitiesArchetype[cur->archetype->entitiesId[i]].index--;
    }
    
    Entity::entitiesArchetype.erase(id);
}

bool Entity::HasComponent(ComponentType compType)
{
    auto archetypeType = Entity::entitiesArchetype[id].archetype->type;
    return std::find(archetypeType.begin(), archetypeType.end(), compType) != archetypeType.end();
}

bool Entity::hasTag(std::string tag)
{
    if (!HasComponent(TagComponent::GetCompId()))
        return false;

    auto tags = GetComponent<TagComponent>()->tags;
    for (size_t i = 0; i < tags.size(); i++)
    {
        if (tags[i] == tag)
            return true;
    }
    return false;
}

void Entity::giveTag(std::string tag) 
{
    TagComponent* tagComp = GetComponent<TagComponent>(true);
    if (tagComp == nullptr)
    {
        tagComp = new TagComponent(id);
        // components.push_back(tagComp);
        AddComponent(tagComp);
    }
    tagComp->tags.push_back(tag);
    if (TagComponent::tagOwnersMap.find(tag) == TagComponent::tagOwnersMap.end())
    {
        TagComponent::tagOwnersMap.emplace(tag, std::unordered_set<unsigned>{id});
    }
    else
    {
        TagComponent::tagOwnersMap.find(tag)->second.emplace(id);
    }
}

Component::Component(unsigned entityId, unsigned componentTypeId)
{
    Component::componentOwnersArray[componentTypeId].emplace(entityId);
    type = componentTypeId;
}

Entity* Entity::EntityById(unsigned Entity_ID, bool is_optional)
{
    auto iterator = Entity::entities_map.find(Entity_ID);
    Entity* res = (iterator == Entity::entities_map.end()) ? nullptr : iterator->second;
    if (res != nullptr || is_optional)
    {
        return res;
    }
    else
    {
        std::cout<<"----------------WRONG ID-----------------\n";
        throw std::exception{};
    }  
}

Entity* Entity::AnyObjectWithComponent(unsigned comp_type,bool can_be_no_objects)
{
    auto iterator = Component::componentOwnersArray[comp_type].begin();
    while (iterator != Component::componentOwnersArray[comp_type].end())
    {
        Entity* entity = Entity::EntityById(*iterator, true);
        if (entity == nullptr)
        {
            iterator = Component::componentOwnersArray[comp_type].erase(iterator);
        }
        else
        {
            return entity;
        }
    }
    if (can_be_no_objects)
    {
        return nullptr;
    }
    else
    {
        std::cout<<"----------------NO SUCH ENTITIES-----------------\n";
        throw std::exception{};
    }
}

Entity* Entity::AnyObjectWithTag(std::string tag,bool optional)
{
    auto iterator_tags = TagComponent::tagOwnersMap.find(tag); 
    if (iterator_tags == TagComponent::tagOwnersMap.end())
    {
        if (optional)
        {
            return nullptr;
        }
        else
        {
            std::cout<<"----------------NO ENTITIES WITH THIS TAG-----------------\n";
            throw std::exception{};
        }
    }
    else
    {
        auto iterator = iterator_tags->second.begin();
        while (iterator != iterator_tags->second.end())
        {
            Entity* entity = Entity::EntityById(*iterator, true);
            if (entity == nullptr)
            {
                iterator = iterator_tags->second.erase(iterator);
            }
            else
            {
                return entity;
            }
        }
        TagComponent::tagOwnersMap.erase(iterator_tags);
        if (optional)
        {
            return nullptr;
        }
        else
        {
            std::cout<<"----------------NO ENTITIES WITH THIS TAG-----------------\n";
            throw std::exception{};
        }
    }
}

std::vector<Entity*> Entity::AllEntitiesWithComponentsAndTags(unsigned n_components, 
    unsigned n_tags, const std::vector<unsigned>& components_v, const std::vector<std::string>& tags_v,
    bool optional)
{
    std::unordered_set<unsigned> result;
    if (n_components > 0)
    {
        result = Component::componentOwnersArray[components_v[0]];
        for (int i = 1; i < n_components; i++)
        {
            intersect(result, Component::componentOwnersArray[components_v[i]]);
        }  
        for (int i = 0; i < n_tags; i++)
        {
            intersect(result, TagComponent::tagOwnersMap[tags_v[i]]);
        }      
    }
    else
    {
        result = TagComponent::tagOwnersMap[tags_v[0]];
        for (int i = 1; i < n_tags; i++)
        {
            intersect(result, TagComponent::tagOwnersMap[tags_v[i]]);
        } 
    }

    std::vector<Entity*> v{};
    for (auto i = result.begin(); i != result.end(); i++)
    {
        Entity* ent = Entity::EntityById(*i, true);
        if (ent == nullptr)
        {
            for (int j = 0; j < n_tags; j++)
            {
                TagComponent::tagOwnersMap[tags_v[j]].erase(*i);
            } 
            for (int j = 0; j < n_components; j++)
            {
                Component::componentOwnersArray[components_v[j]].erase(*i);
            }   
        }
        else
        {
            v.push_back(ent);
        }
    }

    if (v.size() == 0)
    {
        if (optional)
        {
            return std::vector<Entity*>{};
        }
        else
        {
            std::cout<<"----------------NO ENTITIES WITH THIS COMPONENTS AND TAGS [AllEntitiesWithComponentsAndTags]-----------------\n";
            throw std::exception{};
        }
    }
    else
    {
        return v;
    }
}

void Entity::AddComponent(Component *comp)
{
    // debug("RIGHT IN AddComponent");
    // if (comp.type == 11)
    // {
    //     debug(static_cast<SunLightComponent*>(&comp)->SHADOW_WIDTH);
    // }
    entitiesArchetype[id].archetype->AddComponentToEntity(id, comp);
}

// std::vector<Component*> Component::AnyComponentsInSameEntityWithTags(unsigned n_components, 
//     unsigned n_tags, const std::vector<unsigned>& components_v, const std::vector<std::string>& tags_v,
//     bool optional)
// {
//     std::unordered_set<unsigned> result;
//     if (n_components > 0)
//     {
//         result = Component::componentOwnersArray[components_v[0]];
//         for (int i = 1; i < n_components; i++)
//         {
//             intersect(result, Component::componentOwnersArray[components_v[i]]);
//         }  
//         for (int i = 0; i < n_tags; i++)
//         {
//             intersect(result, TagComponent::tagOwnersMap[tags_v[i]]);
//         }      
//     }
//     else
//     {
//         result = TagComponent::tagOwnersMap[tags_v[0]];
//         for (int i = 1; i < n_tags; i++)
//         {
//             intersect(result, TagComponent::tagOwnersMap[tags_v[i]]);
//         } 
//     }

//     std::vector<Component*> v(n_components);
//     for (auto i = result.begin(); i != result.end(); i++)
//     {
//         Entity* ent = Entity::EntityById(*i, true);
//         if (ent == nullptr)
//         {
//             for (int j = 0; j < n_tags; j++)
//             {
//                 TagComponent::tagOwnersMap[tags_v[j]].erase(*i);
//             } 
//             for (int j = 0; j < n_components; j++)
//             {
//                 Component::componentOwnersArray[components_v[j]].erase(*i);
//             }   
//         }
//         else
//         {
//             std::vector<unsigned> compPositions(n_components);
//             ent->HasComponent_OLD(n_components, components_v.data(), compPositions.data());
//             for (int j = 0; j < compPositions.size(); j++)
//             {
//                 v[j] = ent->components[compPositions[j]];
//             }
//             return v;
//         }
//     }
//     if (optional)
//     {
//         return std::vector<Component*>{};
//     }
//     else
//     {
//         std::cout<<"----------------NO ENTITIES WITH THIS COMPONENTS AND TAGS [AnyComponentsInSameEntityWithTags]-----------------\n";
//         throw std::exception{};
//     }
// }

Archetype::Archetype(std::vector<ComponentType> componentList)
{
    type = componentList;
    for (int i = 0; i < type.size(); i++)
    {
        components.push_back(std::vector<Component*>());
    }
    hasUnhandledArchetypes = true;
}

int Archetype::AddEmptyEntity(unsigned eid)
{
    if (!basicNode)
    {
        basicNode = new Archetype(std::vector<ComponentType>{});
    }
    basicNode->entitiesId.push_back(eid);
    Entity::entitiesArchetype.emplace(eid, Record(basicNode, (unsigned)(basicNode->entitiesId.size() - 1)));
    return basicNode->entitiesId.size();
}

void Archetype::AddComponentToEntity(unsigned int EntityIndex, Component *component)
{
    if (IS_ARCHETYPE_DEBUG)
        std::cout << "ADDING COMPONENT: " << component->type << std::endl;
    int i = 0;
    for (; i < mates.size(); i++)
    {
        if(mates[i].componentLink == component->type)
            break;
    }
    if (i == mates.size())
    {
        const auto recursiveFinder = [](std::vector<ComponentType> components, 
                                        Archetype*  startingNode,
                                        const auto& finderFunction) -> Archetype*
        {
            if (IS_ARCHETYPE_DEBUG)
            {
                startingNode->debugPrintType();
                debug(components.size());
            }

            if (components.size() == 0)
                return startingNode;

            for (int i = 0; i < startingNode->mates.size(); i++)
            {
                auto iter = std::find(components.begin(), components.end(), startingNode->mates[i].componentLink);
                if (iter != components.end())
                {
                    std::vector<ComponentType> componentsCopy = components;
                    Archetype *newStartingNode = startingNode->mates[i].mate;
                    componentsCopy.erase(iter - components.begin() + componentsCopy.begin());
                    Archetype* result = 
                        finderFunction(std::move(componentsCopy), newStartingNode, finderFunction);
                    if (result != nullptr)
                    {
                        if (IS_ARCHETYPE_DEBUG)
                            debug("success");
                        return result;
                    }
                }
            }
            if (IS_ARCHETYPE_DEBUG)
                debug("null");
            return nullptr;  
        };
        std::vector<ComponentType> componentList = type;
        componentList.push_back(component->type);
        Archetype* result = recursiveFinder(componentList, Archetype::basicNode, recursiveFinder);
        if (result == nullptr)
        {
            result = new Archetype(componentList);
        }
        ArchetypeGraphEdge edge;
        edge.mate = result;
        edge.componentLink = component->type;
        mates.push_back(std::move(edge));
    }
    Archetype* newArchetype = mates[i].mate;
    Record& record = Entity::entitiesArchetype[EntityIndex];
    std::vector<std::vector<Component*>>* currentComponentsData = &(record.archetype->components);
    for (int j = 0; j < currentComponentsData->size(); j++)
    {
        for (int k = 0; k < newArchetype->type.size(); k++)
        {
            if(newArchetype->type[k] == (*currentComponentsData)[j][0]->type)
            {
                newArchetype->components[k].push_back((*currentComponentsData)[j][record.index]);
                break;
            }
        }
        (*currentComponentsData)[j].erase((*currentComponentsData)[j].begin() + record.index);
    }
    for (int k = 0; k < newArchetype->type.size(); k++)
    {
        if(newArchetype->type[k] == component->type)
        {
            newArchetype->components[k].push_back(component);
            break;
        }
    }
    newArchetype->entitiesId.push_back(entitiesId[record.index]);
    entitiesId.erase(entitiesId.begin() + record.index);
    for (int i = record.index; i < record.archetype->entitiesId.size(); i++)
    {
        Entity::entitiesArchetype[record.archetype->entitiesId[i]].index--;
    }
    record.archetype = newArchetype;
    record.index = newArchetype->components[0].size() - 1;
}

void Archetype::debugPrintType()
{
    std::cout << "TYPE: ";
    for (int i = 0; i < type.size(); i++)
    {
        std::cout << type[i] << " ";
    }
    std::cout << "\n";
}

void Archetype::debugShowProducedTree()
{
    if (this == Archetype::basicNode)
        std::cout << "ROOT\n";
    else
        std::cout << "NODE\n";
    debugPrintType();
    std::cout << "Entities: " << entitiesId.size() << std::endl;
    std::cout << "Connections: ";
    for (int i = 0; i < mates.size(); i++)
    {
        std::cout << mates[i].componentLink << " ";          
    }
    std::cout << std::endl;
    for (int i = 0; i < mates.size(); i++)
    {
        mates[i].mate->debugShowProducedTree();         
    }
}

void Archetype::debugShowTree()
{
    if (basicNode == nullptr)
    {
        debug("EMPTY TREE");
        return;
    }
    Archetype::basicNode->debugShowProducedTree();
}

System::System()
{
    allSystems.push_back(this);
    this->UpdateArchetypes();
    isSeparateSystem = true;
}

void System::UpdateArchetypes()
{
    if (activeComponents.size() == 0)
        return;
    

    const auto recursiveAppender = [this](std::vector<ComponentType> components, 
                                        Archetype*  startingNode,
                                        const auto& appenderFunction) -> void
    {
        if (!startingNode)
            return;

        if (components.size() == 0)
        {
            relevantArchetypes.emplace(startingNode);
            for (int i = 0; i < startingNode->mates.size(); i++)
            {
                appenderFunction(std::vector<ComponentType>{},
                                 startingNode->mates[i].mate,
                                 appenderFunction);
            }
        }
        else
        {
            for (int i = 0; i < startingNode->mates.size(); i++)
            {
                std::vector<ComponentType> componentsCopy = components;
                Archetype *newStartingNode = startingNode->mates[i].mate;
                auto iter = std::find(components.begin(), components.end(), startingNode->mates[i].componentLink);
                if (iter != components.end())
                    componentsCopy.erase(iter - components.begin() + componentsCopy.begin());
                appenderFunction(std::move(componentsCopy), newStartingNode, appenderFunction);
            }
        }
    };
    recursiveAppender(activeComponents, Archetype::basicNode, recursiveAppender);    
}

void System::UpdateAllSystems()
{
    for (int i = 0; i < allSystems.size(); i++)
    {
        allSystems[i]->UpdateArchetypes();
    }   
}

void Archetype::handleNewArchetypes()
{
    if (hasUnhandledArchetypes)
    {
        System::UpdateAllSystems();
        hasUnhandledArchetypes = false;
    }
}

void System::RunAllSystems()
{
    for (int i = 0; i < allSystems.size(); i++)
    {
        System* system = allSystems[i];
        
        if (!system->isSeparateSystem)
            continue;

        for (auto &it : system->relevantArchetypes)
        {
            std::vector<unsigned> positionsInArchetype{};
            for (int j = 0; j < system->activeComponents.size(); j++)
            {
                ComponentType cur = system->activeComponents[j];
                unsigned pos = std::find(it->type.begin(),it->type.end(), cur) - it->type.begin();
                positionsInArchetype.push_back(pos);
            }
            for (int j = 0; j < it->entitiesId.size(); j++)
            {
                std::vector<Component*> systemComponents;
                for (int k = 0; k < positionsInArchetype.size(); k++)
                {
                    systemComponents.push_back(it->components[positionsInArchetype[k]][j]);
                }
                system->RunSystem(systemComponents);
            }
        }    
    } 
}

void System::AddQueue(std::string name, std::vector<ComponentType> comps)
{
    Queue* queue = new Queue(name);
    queue->activeComponents = comps;
    queues.emplace(queue);
    queue->UpdateArchetypes();
}

Queue* System::GetQueueByName(std::string name)
{
    for (auto &it: queues)
    {
        if (it->name == name)
            return it;
    }
    debug("NO SUCH QUEUE!!!");
    throw std::exception{};
    return nullptr;
}

void Queue::RunSystem(std::vector<Component*> comp)
{
    throw std::exception{};
}

Queue::Queue(std::string s) : name{s}
{
    isSeparateSystem = false;
}

Queue::Iterator::const_reference Queue::Iterator::operator*() const
{
    return m_currentComponents;
}

bool Queue::Iterator::IsEndIterator()
{
    return m_archetypesIter == m_parentQueue->relevantArchetypes.end();
}

void Queue::Iterator::UpdatePositionsInArchetype()
{
    Archetype* archetype = *m_archetypesIter;
    m_componentPositionsInArchetype.clear();
    for (int j = 0; j < m_parentQueue->activeComponents.size(); j++)
    {
        ComponentType cur = m_parentQueue->activeComponents[j];
        unsigned pos = std::find(archetype->type.begin(),archetype->type.end(), cur) 
                            - archetype->type.begin();
        m_componentPositionsInArchetype.push_back(pos);
    }
}

bool Queue::Iterator::ValidateIndexes()
{
    bool flag = false;
    while (m_archetypesIter != m_parentQueue->relevantArchetypes.end())
    {
        Archetype* archetype = *m_archetypesIter;
        if (m_entityInArchetypeNum >= archetype->entitiesId.size())
        {
            m_entityInArchetypeNum -= archetype->entitiesId.size();
            m_archetypesIter++;
            flag = true;
        }
        else
        {
            break;
        }
    }
    return flag;
}

void Queue::Iterator::UpdateCurrentComponents()
{
    if (IsEndIterator())
        return;

    m_currentComponents.clear();
    Archetype* archetype = *m_archetypesIter;
    for (int k = 0; k < m_componentPositionsInArchetype.size(); k++)
    {
        m_currentComponents.push_back(
            (*m_archetypesIter)->components[m_componentPositionsInArchetype[k]][m_entityInArchetypeNum]);
    }
}

Queue::Iterator& Queue::Iterator::operator++()
{
    m_entityInArchetypeNum++;
    if (ValidateIndexes())
    {
        if (!IsEndIterator())
            UpdatePositionsInArchetype();
    }
    UpdateCurrentComponents();
    return *this;
}

Queue::Iterator Queue::Iterator::operator++(int)
{
    auto tmp = *this;
    ++(*this);
    return tmp;
}

bool operator== (const Queue::Iterator& a, const Queue::Iterator& b)
{
    if (a.m_parentQueue != b.m_parentQueue)
    {
        return false;
    }
    if (a.m_archetypesIter != b.m_archetypesIter)
    {
        return false;
    }
    if (a.m_archetypesIter == a.m_parentQueue->relevantArchetypes.end())
    {
        return true;
    }
    return a.m_entityInArchetypeNum == b.m_entityInArchetypeNum;
}

bool operator!= (const Queue::Iterator& a, const Queue::Iterator& b)
{
    return !(a == b);
}

Queue::Iterator::Iterator(Queue* queue, bool is_end = false)
{
    m_currentComponents = value_type{};
    m_entityInArchetypeNum = 0;
    m_componentPositionsInArchetype = std::vector<unsigned> {};
    m_parentQueue = queue;

    if (!is_end)
    {
        m_archetypesIter = queue->relevantArchetypes.begin();
        ValidateIndexes();
        if (!IsEndIterator())
            UpdatePositionsInArchetype();
        UpdateCurrentComponents();
    }
    else
    {
        m_archetypesIter = queue->relevantArchetypes.end();
    }
}

Queue::Iterator Queue::begin()
{
    return Iterator(this);
}

Queue::Iterator Queue::end()
{
    return Iterator(this, true);
}